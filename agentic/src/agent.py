"""
GeneRNBI orchestrating agent — powered by Biomni A1.

Two tools:
  - search_docs:        RAG over agentic/docs/ (API spec, datasets, integration guide)
  - search_manuscript:  RAG over the benchmark manuscript PDF (paper content)
"""

import os
from pathlib import Path
from typing import Optional

from dotenv import load_dotenv

from config import LLM_PROVIDER, MODEL_ID
from rag_builder import (
    _missing_key_error,
    get_or_build_docs_index,
    get_or_build_manuscript_index,
)

_SOURCE_MAP = {"openai": "OpenAI", "anthropic": "Anthropic"}


# ─── Tool schemas (Biomni format) ────────────────────────────────────────────

_SEARCH_DOCS_SCHEMA = {
    "name": "search_docs",
    "description": (
        "Search the GeneRNBI documentation for information about integrating new methods, "
        "metrics, or datasets; the API format for Viash components (config.vsh.yaml, "
        "script.py structure, argument types, file formats); dataset descriptions; "
        "and any other framework usage or integration question. "
        "Use this tool FIRST for all technical and integration questions."
    ),
    "required_parameters": [
        {
            "name": "query",
            "type": "str",
            "description": "Any question about framework usage, integration, or API structure.",
            "default": None,
        }
    ],
    "optional_parameters": [],
}

_SEARCH_MANUSCRIPT_SCHEMA = {
    "name": "search_manuscript",
    "description": (
        "Search the GeneRNBI benchmark manuscript PDF for paper-level content: "
        "experimental results, biological motivation, method comparisons, design "
        "decisions, and conclusions. "
        "Use this tool ONLY when search_docs cannot answer — i.e. when the question "
        "explicitly concerns the paper's findings or scientific rationale."
    ),
    "required_parameters": [
        {
            "name": "query",
            "type": "str",
            "description": "Question about the manuscript content.",
            "default": None,
        }
    ],
    "optional_parameters": [],
}


# ─── Tool registration helper ─────────────────────────────────────────────────

def _register_tool(a1, fn, schema: dict):
    """Directly register a tool with A1, bypassing LLM-based schema generation."""
    import builtins

    name = schema["name"]
    module = "grn_tools"

    if not hasattr(a1, "module2api") or a1.module2api is None:
        a1.module2api = {}
    a1.module2api.setdefault(module, [])
    a1.module2api[module] = [t for t in a1.module2api[module] if t.get("name") != name]
    a1.module2api[module].append({**schema, "module": module})

    if not hasattr(a1, "_custom_functions"):
        a1._custom_functions = {}
    a1._custom_functions[name] = fn

    if not hasattr(a1, "_custom_tools"):
        a1._custom_tools = {}
    a1._custom_tools[name] = {"name": name, "description": schema["description"], "module": module}

    if not hasattr(builtins, "_biomni_custom_functions"):
        builtins._biomni_custom_functions = {}
    builtins._biomni_custom_functions[name] = fn

    print(f"Registered tool '{name}'")


# ─── Tool implementations ─────────────────────────────────────────────────────

def _build_docs_search_fn(index, li_llm):
    """Direct retriever for docs — returns raw chunks without LLM synthesis.
    
    Since docs are small (2 files, ~430 lines), returning the top chunks verbatim
    gives more accurate, detailed answers than LLM synthesis which tends to summarize.
    """
    retriever = index.as_retriever(similarity_top_k=10)

    def _search(query: str) -> str:
        try:
            nodes = retriever.retrieve(query)
            if not nodes:
                result = "No relevant content found in docs."
            else:
                parts = [node.get_content() for node in nodes]
                result = "\n\n---\n\n".join(parts)
            print(result)
            return result
        except Exception as e:
            err = f"Error searching docs: {e}"
            print(err)
            return err

    return _search


def _build_manuscript_search_fn(index, li_llm, n_subqueries: int = 3):
    """Parallel multi-query search for manuscript (large PDF, sub-query expansion helps)."""
    from concurrent.futures import ThreadPoolExecutor, as_completed
    from llama_index.core.prompts import PromptTemplate

    query_engine = index.as_query_engine(
        similarity_top_k=8,
        llm=li_llm,
        response_mode="compact",
    )

    _EXPAND_TMPL = PromptTemplate(
        f"Break the following question into exactly {n_subqueries} precise, "
        "non-overlapping sub-questions that together give a complete answer.\n"
        "Return ONLY the sub-questions, one per line, no numbering or prefixes.\n\n"
        "Question: {question}\n\nSub-questions:"
    )

    def _search(query: str) -> str:
        try:
            expansion = li_llm.complete(_EXPAND_TMPL.format(question=query)).text
            sub_questions = [q.strip() for q in expansion.strip().splitlines() if q.strip()][:n_subqueries]

            def _query(q):
                return (q, str(query_engine.query(q)))

            parts: list = []
            seen: set = set()
            with ThreadPoolExecutor(max_workers=n_subqueries) as pool:
                futures = {pool.submit(_query, q): q for q in sub_questions}
                for future in as_completed(futures):
                    q, resp = future.result()
                    key = resp[:120]
                    if key not in seen:
                        seen.add(key)
                        parts.append(f"[{q}]\n{resp}")

            result = "\n\n---\n\n".join(parts) if parts else "No relevant content found."
            print(result)
            return result
        except Exception as e:
            err = f"Error searching manuscript: {e}"
            print(err)
            return err

    return _search


# ─── System prompt ────────────────────────────────────────────────────────────

_SYSTEM_PROMPT_PREFIX = """\
You are the GeneRNBI Assistant — the primary interface for the GeneRNBI GRN benchmark framework.

══════════════════════════════════════════════════════
STRICT RULES — FOLLOW THESE EXACTLY
══════════════════════════════════════════════════════

1. You have EXACTLY TWO tools: search_docs and search_manuscript.
   - ALWAYS call search_docs FIRST for ANY question about datasets, methods,
     metrics, file formats, or integration.
   - ONLY call search_manuscript when the user asks about the paper's findings,
     results, or scientific motivation.

2. NEVER use os, sys, glob, subprocess, pathlib, or any Python library to answer.
   Do NOT list files, read files, or access the filesystem.

3. NEVER answer from memory. ALWAYS call a tool first.

4. For EVERY question, your FIRST action must be:
     result = search_docs('...your query...'); print(result)
   Then synthesise the answer from the tool output.

══════════════════════════════════════════════════════
TOOL USAGE
══════════════════════════════════════════════════════

search_docs  — Use for ALL technical questions:
  • Datasets available (names, formats, which metrics apply)
  • How to add a new method or metric (file structure, config.vsh.yaml format)
  • Viash component format (arguments, types, directions, __merge__ field)
  • API spec and file formats (h5ad, csv, gmt, etc.)
  • Which GRN methods and metrics are currently available
  Usage:  result = search_docs('your question'); print(result)

search_manuscript  — ONLY for paper-level content:
  • Experimental results and benchmarking findings
  • Biological motivation and scientific rationale
  • Performance comparisons between methods
  • Design decisions and conclusions in the paper
  Usage:  result = search_manuscript('your question'); print(result)

"""


# ─── Main initializer ─────────────────────────────────────────────────────────

def initialize_agent(
    docs_dir: Optional[Path] = None,
    data_dir: Optional[Path] = None,
    use_cache: bool = True,
    model_id: str = None,
    llm_provider: str = None,
):
    """Build indices and return an A1 orchestrator with two RAG tools."""
    load_dotenv()

    provider = (llm_provider or LLM_PROVIDER).lower()
    model = model_id or MODEL_ID
    source = _SOURCE_MAP.get(provider)
    if source is None:
        raise ValueError(f"Unsupported provider '{provider}'. Choose 'openai' or 'anthropic'.")

    if provider == "openai" and not os.environ.get("OPENAI_API_KEY"):
        _missing_key_error()
    elif provider == "anthropic" and not os.environ.get("ANTHROPIC_API_KEY"):
        _missing_key_error()

    agentic_root = Path(__file__).parent.parent
    if docs_dir is None:
        docs_dir = agentic_root / "docs"
    if data_dir is None:
        data_dir = agentic_root / "data"

    print("Building docs index...")
    docs_index = get_or_build_docs_index(docs_dir, use_cache=use_cache)

    print("Building manuscript index...")
    manuscript_index = get_or_build_manuscript_index(data_dir, use_cache=use_cache)

    # LlamaIndex LLM for query engines and sub-query expansion
    if provider == "openai":
        from llama_index.llms.openai import OpenAI as LIOpenAI
        _li_llm = LIOpenAI(model=model)
    else:
        from llama_index.llms.anthropic import Anthropic as LIAnthropic
        _li_llm = LIAnthropic(model=model)

    print("Creating A1 orchestrator...")
    from biomni.agent.a1 import A1

    a1 = A1(
        path=str(agentic_root / ".biomni_data"),
        llm=model,
        source=source,
        use_tool_retriever=False,
        expected_data_lake_files=[],
    )

    print("Registering tools...")
    a1.module2api = {}

    search_docs_fn = _build_docs_search_fn(docs_index, _li_llm)
    search_docs_fn.__doc__ = (
        "Search GeneRNBI documentation for integration guides, API formats, "
        "Viash component structure, dataset info, and framework usage.\n\n"
        "Args:\n    query (str): Any technical or integration question.\n\n"
        "Returns:\n    str: Answer from the documentation."
    )
    search_docs_fn.__name__ = "search_docs"

    search_manuscript_fn = _build_manuscript_search_fn(manuscript_index, _li_llm)
    search_manuscript_fn.__doc__ = (
        "Search the GeneRNBI benchmark manuscript for paper-level content: "
        "results, motivation, comparisons, design decisions, conclusions.\n\n"
        "Args:\n    query (str): Question about the manuscript content.\n\n"
        "Returns:\n    str: Answer from the manuscript."
    )
    search_manuscript_fn.__name__ = "search_manuscript"

    _register_tool(a1, search_docs_fn, _SEARCH_DOCS_SCHEMA)
    _register_tool(a1, search_manuscript_fn, _SEARCH_MANUSCRIPT_SCHEMA)

    a1.configure()
    a1._inject_custom_functions_to_repl()

    a1.system_prompt = _SYSTEM_PROMPT_PREFIX + a1.system_prompt

    print("A1 agent ready!")
    return a1


