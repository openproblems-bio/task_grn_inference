"""
GeneRNBI orchestrating agent — powered by Biomni A1.

The A1 agent orchestrates two specialised tools:
  - search_manuscript:       RAG over the benchmark manuscript PDF
  - call_integration_agent:  LlamaIndex sub-agent with code/docs/viash access
"""

import asyncio
import os
from pathlib import Path
from typing import Optional

from dotenv import load_dotenv

from config import LLM_PROVIDER, MODEL_ID
from rag_builder import _missing_key_error, get_or_build_docs_index, get_or_build_manuscript_index
from integration_agent import create_integration_agent

# Map our config provider names to biomni SourceType values
_SOURCE_MAP = {"openai": "OpenAI", "anthropic": "Anthropic"}

_SEARCH_MANUSCRIPT_SCHEMA = {
    "name": "search_manuscript",
    "description": (
        "Search the GeneRNBI benchmark manuscript PDF for information about metrics, "
        "methods, datasets, pipeline design, standardization, experimental results, "
        "or any other paper content. For broad questions the tool internally expands "
        "the query into focused sub-questions and returns a comprehensive answer."
    ),
    "required_parameters": [
        {
            "name": "query",
            "type": "str",
            "description": "Any question about the manuscript — broad or specific.",
            "default": None,
        }
    ],
    "optional_parameters": [],
}

_INTEGRATION_AGENT_SCHEMA = {
    "name": "call_integration_agent",
    "description": (
        "Call the GeneRNBI Integration Agent for code-related tasks: integrating new GRN "
        "inference methods (src/methods/), adding evaluation metrics (src/metrics/), "
        "adding datasets, reading or debugging source code, configuring Viash components "
        "(config.vsh.yaml), and troubleshooting pipeline errors or layer mismatches. "
        "The agent has direct access to source code and live Viash documentation at viash.io."
    ),
    "required_parameters": [
        {
            "name": "query",
            "type": "str",
            "description": "Integration or troubleshooting question.",
            "default": None,
        }
    ],
    "optional_parameters": [],
}


def _register_tool(a1, fn, schema: dict):
    """Directly register a tool with A1, bypassing the LLM-based schema generation."""
    import builtins

    name = schema["name"]
    module = "grn_tools"

    if not hasattr(a1, "module2api") or a1.module2api is None:
        a1.module2api = {}
    a1.module2api.setdefault(module, [])
    # Remove any existing entry with same name
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



def _build_search_manuscript_fn(manuscript_index, llm_complete_fn):
    """Return a sync search_manuscript callable bound to the given index."""
    from llama_index.core.response_synthesizers import get_response_synthesizer
    from llama_index.core.prompts import PromptTemplate

    base_engine = manuscript_index.as_query_engine(
        similarity_top_k=15,
        response_synthesizer=get_response_synthesizer(response_mode="tree_summarize"),
    )

    _EXPAND_TMPL = PromptTemplate(
        "You are a research assistant. Break the following question into 4-6 precise, "
        "non-overlapping sub-questions that together give a COMPLETE answer.\n"
        "Return ONLY the sub-questions, one per line, no numbering.\n\n"
        "Question: {question}\n\nSub-questions:"
    )

    def search_manuscript(query: str) -> str:
        """Search the GeneRNBI benchmark manuscript PDF for information about metrics,
        methods, datasets, pipeline design, standardization, experimental results, and
        any other content from the paper.

        For broad questions the tool internally expands the query into focused
        sub-questions and returns a comprehensive synthesised answer.

        Args:
            query (str): Any question about the manuscript — broad or specific.

        Returns:
            str: Comprehensive answer drawn from the manuscript.
        """
        try:
            expansion = llm_complete_fn(_EXPAND_TMPL.format(question=query))
            sub_questions = [q.strip() for q in expansion.strip().splitlines() if q.strip()]
            queries = [query] + sub_questions[:5]

            parts: list = []
            seen: set = set()
            for q in queries:
                resp = str(base_engine.query(q))
                key = resp[:120]
                if key not in seen:
                    seen.add(key)
                    parts.append(f"[Sub-query: {q}]\n{resp}")

            return "\n\n---\n\n".join(parts)
        except Exception as e:
            return f"Error searching manuscript: {e}"

    return search_manuscript


def _build_call_integration_agent_fn(integration_ag):
    """Return a sync call_integration_agent callable wrapping the async sub-agent."""

    def call_integration_agent(query: str) -> str:
        """Call the GeneRNBI Integration Agent for code-related tasks including:
        integrating new GRN inference methods (src/methods/), adding evaluation
        metrics (src/metrics/), adding datasets, reading or debugging source code,
        configuring Viash components (config.vsh.yaml), and troubleshooting
        pipeline errors or layer mismatches.

        The integration agent has direct access to the source code and live
        Viash documentation at viash.io.

        Args:
            query (str): Integration or troubleshooting question.

        Returns:
            str: Answer from the integration agent.
        """
        try:
            async def _run():
                handler = integration_ag.run(query)
                return await handler

            return str(asyncio.run(_run()))
        except RuntimeError:
            # Already inside an event loop — use a new thread
            import concurrent.futures
            with concurrent.futures.ThreadPoolExecutor(max_workers=1) as pool:
                future = pool.submit(asyncio.run, _run())
                return str(future.result())
        except Exception as e:
            return f"Integration agent error: {e}"

    return call_integration_agent


def initialize_agent(
    docs_dir: Optional[Path] = None,
    data_dir: Optional[Path] = None,
    use_cache: bool = True,
    model_id: str = None,
    llm_provider: str = None,
):
    """Build all indices, create sub-agents, and return an A1 orchestrator."""
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

    print("Creating integration agent...")
    integration_ag = create_integration_agent(docs_index, model_id=model, llm_provider=provider)

    # Build A1 orchestrator
    print("Creating A1 orchestrator...")
    from biomni.agent.a1 import A1

    a1 = A1(
        path=str(agentic_root / ".biomni_data"),
        llm=model,
        source=source,
        use_tool_retriever=False,
        expected_data_lake_files=[],  # skip biomni data lake download
    )

    # Build a lightweight LLM complete function for query expansion
    if provider == "openai":
        from llama_index.llms.openai import OpenAI as LIOpenAI
        _li_llm = LIOpenAI(model=model)
    else:
        from llama_index.llms.anthropic import Anthropic as LIAnthropic
        _li_llm = LIAnthropic(model=model)

    def _llm_complete(prompt_str: str) -> str:
        return _li_llm.complete(prompt_str).text

    # Register tools — clear default biomni tools first to keep system prompt compact
    print("Registering tools...")
    a1.module2api = {}  # remove all biomni built-in tools
    search_fn = _build_search_manuscript_fn(manuscript_index, _llm_complete)
    integration_fn = _build_call_integration_agent_fn(integration_ag)

    _register_tool(a1, search_fn, _SEARCH_MANUSCRIPT_SCHEMA)
    _register_tool(a1, integration_fn, _INTEGRATION_AGENT_SCHEMA)
    a1.configure()  # rebuild LangGraph with updated module2api
    a1._inject_custom_functions_to_repl()  # make fns callable in Python REPL

    # Prepend GRN-specific routing instructions to the generated system prompt
    instructions = (Path(__file__).parent / "instructions.txt").read_text()
    no_import_note = (
        "IMPORTANT: The functions `search_manuscript` and `call_integration_agent` are "
        "pre-loaded in the execution environment. Do NOT use any import statement for them. "
        "Call them directly, e.g.:\n"
        "  result = search_manuscript('your question here'); print(result)\n"
        "  result = call_integration_agent('your question here'); print(result)\n\n"
    )
    a1.system_prompt = instructions + "\n\n" + no_import_note + a1.system_prompt

    print("A1 agent ready!")
    return a1
