"""
GeneRNBI orchestrating agent — powered by Biomni A1.

The A1 agent orchestrates three specialised tools:
  - read_file:          read any source file in task_grn_inference/
  - list_directory:     list files/subdirs in any allowed path
  - search_manuscript:  parallel RAG over the benchmark manuscript PDF
"""

import os
from pathlib import Path
from typing import Optional

from dotenv import load_dotenv

from config import LLM_PROVIDER, MODEL_ID
from rag_builder import _missing_key_error, get_or_build_manuscript_index

# Map our config provider names to biomni SourceType values
_SOURCE_MAP = {"openai": "OpenAI", "anthropic": "Anthropic"}

# Root of the task_grn_inference repository (one level above agentic/)
_REPO_ROOT = Path(__file__).parent.parent.parent.resolve()


# ─── Tool schemas (Biomni format) ────────────────────────────────────────────

_READ_FILE_SCHEMA = {
    "name": "read_file",
    "description": (
        "Read the full contents of a source file inside the task_grn_inference "
        "repository. Use this to inspect scripts (script.py), Viash component "
        "configs (config.vsh.yaml), utility modules, metric implementations, or "
        "inference method code. Paths are relative to the task_grn_inference root, "
        "e.g. 'src/metrics/tf_binding/script.py'."
    ),
    "required_parameters": [
        {
            "name": "path",
            "type": "str",
            "description": "File path relative to task_grn_inference root.",
            "default": None,
        }
    ],
    "optional_parameters": [],
}

_LIST_DIRECTORY_SCHEMA = {
    "name": "list_directory",
    "description": (
        "List files and subdirectories inside the task_grn_inference repository. "
        "Use this to discover what methods (src/methods/), metrics (src/metrics/), "
        "or other components exist. Paths are relative to the task_grn_inference root. "
        "Pass '' or '.' to list the repo root."
    ),
    "required_parameters": [
        {
            "name": "path",
            "type": "str",
            "description": "Directory path relative to task_grn_inference root.",
            "default": None,
        }
    ],
    "optional_parameters": [],
}

_SEARCH_MANUSCRIPT_SCHEMA = {
    "name": "search_manuscript",
    "description": (
        "Search the GeneRNBI benchmark manuscript PDF for information not present "
        "in the source code: experimental results, biological motivation, comparisons "
        "between methods, design decisions, and paper-level conclusions. "
        "Use read_file / list_directory FIRST for any code or config question. "
        "Only call this tool when the answer requires paper content."
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

def _build_read_file_fn(repo_root: Path):
    def read_file(path: str) -> str:
        """Read the contents of a source file in task_grn_inference.

        Args:
            path (str): File path relative to the task_grn_inference root,
                        e.g. 'src/metrics/tf_binding/script.py'.

        Returns:
            str: File contents, or an error message if the file is not found.
        """
        try:
            target = (repo_root / path).resolve()
            if not str(target).startswith(str(repo_root)):
                result = "Error: path is outside the task_grn_inference repository."
            elif not target.exists():
                result = f"File not found: {path}"
            elif not target.is_file():
                result = f"Not a file: {path} — use list_directory to explore directories."
            else:
                result = target.read_text(encoding="utf-8", errors="replace")
        except Exception as e:
            result = f"Error reading file: {e}"
        print(result)
        return result

    return read_file


def _build_list_directory_fn(repo_root: Path):
    def list_directory(path: str) -> str:
        """List files and subdirectories in task_grn_inference.

        Args:
            path (str): Directory path relative to the task_grn_inference root.
                        Pass '' or '.' to list the repo root.

        Returns:
            str: Newline-separated list of entries, or an error message.
        """
        try:
            target = (repo_root / path).resolve() if path.strip() not in ("", ".") else repo_root
            if not str(target).startswith(str(repo_root)):
                result = "Error: path is outside the task_grn_inference repository."
            elif not target.exists():
                result = f"Directory not found: {path}"
            elif not target.is_dir():
                result = f"Not a directory: {path} — use read_file to read files."
            else:
                entries = sorted(target.iterdir(), key=lambda p: (p.is_file(), p.name))
                lines = [f"{e.name}{'/' if e.is_dir() else ''}" for e in entries]
                result = "\n".join(lines) if lines else "(empty directory)"
        except Exception as e:
            result = f"Error listing directory: {e}"
        print(result)
        return result

    return list_directory


def _build_search_manuscript_fn(manuscript_index, li_llm):
    """Return a fast, parallel search_manuscript callable."""
    from concurrent.futures import ThreadPoolExecutor, as_completed
    from llama_index.core import Settings
    from llama_index.core.prompts import PromptTemplate

    query_engine = manuscript_index.as_query_engine(
        similarity_top_k=8,
        llm=li_llm,
        response_mode="compact",
    )

    _EXPAND_TMPL = PromptTemplate(
        "Break the following question into exactly 3 precise, non-overlapping sub-questions "
        "that together give a complete answer.\n"
        "Return ONLY the sub-questions, one per line, no numbering or prefixes.\n\n"
        "Question: {question}\n\nSub-questions:"
    )

    def search_manuscript(query: str) -> str:
        """Search the GeneRNBI benchmark manuscript PDF.

        Use only for paper-level content: experimental results, biological
        motivation, method comparisons, design decisions, and conclusions.
        For code or config questions, use read_file or list_directory instead.

        Args:
            query (str): Question about the manuscript content.

        Returns:
            str: Answer synthesised from the manuscript.
        """
        try:
            expansion = li_llm.complete(_EXPAND_TMPL.format(question=query)).text
            sub_questions = [q.strip() for q in expansion.strip().splitlines() if q.strip()][:3]

            def _query(q):
                return (q, str(query_engine.query(q)))

            parts: list = []
            seen: set = set()
            with ThreadPoolExecutor(max_workers=3) as pool:
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
            return f"Error searching manuscript: {e}"

    return search_manuscript


# ─── System prompt content ────────────────────────────────────────────────────

def _build_system_prompt_prefix(docs_dir: Path) -> str:
    """Build a rich framework-overview prefix for A1's system prompt."""
    readme = (docs_dir / "README.md").read_text(encoding="utf-8") if (docs_dir / "README.md").exists() else ""
    datasets = (docs_dir / "datasets.md").read_text(encoding="utf-8") if (docs_dir / "datasets.md").exists() else ""

    return f"""You are the GeneRNBI Assistant — the primary interface for the geneRNIB GRN benchmark framework.

══════════════════════════════════════════════════════════════════
FRAMEWORK DOCUMENTATION
══════════════════════════════════════════════════════════════════

{readme}

──────────────────────────────────────────────────────────────────
DATASETS
──────────────────────────────────────────────────────────────────

{datasets}

══════════════════════════════════════════════════════════════════
REPOSITORY LAYOUT  (task_grn_inference/)
══════════════════════════════════════════════════════════════════

src/
  methods/          — one subdirectory per GRN inference method
                      each contains: script.py  config.vsh.yaml  run.sh
  metrics/          — one subdirectory per evaluation metric
                      each contains: script.py  config.vsh.yaml  helper.py  run.sh
  utils/            — shared utility code (util.py, helper.py)
  api/              — Viash API specs (comp_method.yaml, comp_metric.yaml)
  workflows/        — pipeline workflow definitions
  process_data/     — data preprocessing scripts

agentic/            — this assistant

AVAILABLE METHODS  (src/methods/):
  celloracle, dictys, figr, geneformer, granie, grnboost,
  negative_control, pearson_corr, portia, positive_control,
  ppcor, scenic, scenicplus, scglue, scgpt, scprint, spearman_corr

AVAILABLE METRICS  (src/metrics/):
  regression, sem, vc, ws_distance, tf_recovery, tf_binding,
  gs_recovery, replicate_consistency, anchor_regression, all_metrics

══════════════════════════════════════════════════════════════════
VIASH COMPONENT STRUCTURE
══════════════════════════════════════════════════════════════════

Every method and metric is a Viash component with two core files:

1. config.vsh.yaml   — declares inputs/outputs, types, and defaults
   - Top section:  __merge__: ../../api/comp_method.yaml  (or comp_metric.yaml)
   - name / namespace / info (label, summary, description)
   - arguments:  list of --name, type, direction (input/output), required, default

2. script.py         — the component logic
   - Reads par dict (from Viash arg parsing)
   - Produces output .h5ad file

Common argument types: file, string, integer, boolean, double
Common directions: input (default), output

══════════════════════════════════════════════════════════════════
YOUR TOOLS AND ROUTING RULES
══════════════════════════════════════════════════════════════════

You have THREE tools — use them in this priority order:

1. list_directory(path)  — discover what files and subdirectories exist.
   Call this FIRST when you need to find a file or understand the structure.
   READ-ONLY — cannot create directories.

2. read_file(path)       — read the full contents of any source file.
   Use for: script.py, config.vsh.yaml, helper.py, util.py, workflow files.
   READ-ONLY — cannot write or modify files.
   Paths are relative to task_grn_inference root, e.g.:
     'src/metrics/tf_binding/script.py'
     'src/methods/scenic/config.vsh.yaml'

3. search_manuscript(query) — FALLBACK ONLY.
   Use only when the answer requires paper content: experimental results,
   biological motivation, method comparisons, design rationale, conclusions.
   Do NOT use for code or configuration questions — use read_file instead.

ROUTING PRIORITY:
  • Code, configs, pipeline structure  → list_directory + read_file
  • Errors, debugging                  → read_file to inspect the failing script
  • Integrating new components         → read_file a similar existing component,
                                         then explain what files the user must create
                                         (do NOT attempt to create files yourself)
  • Paper results, motivation, design  → search_manuscript

IMPORTANT: The functions read_file, list_directory, and search_manuscript are
pre-loaded in the execution environment. Do NOT import them. Call directly:
  result = read_file('src/metrics/tf_binding/script.py'); print(result)
  result = list_directory('src/metrics'); print(result)
  result = search_manuscript('how are metrics standardized?'); print(result)

"""


# ─── Main initializer ─────────────────────────────────────────────────────────

def initialize_agent(
    docs_dir: Optional[Path] = None,
    data_dir: Optional[Path] = None,
    use_cache: bool = True,
    model_id: str = None,
    llm_provider: str = None,
):
    """Build all indices and return an A1 orchestrator with GRN tools."""
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

    print("Building manuscript index...")
    manuscript_index = get_or_build_manuscript_index(data_dir, use_cache=use_cache)

    # Build LlamaIndex LLM for manuscript query engine
    if provider == "openai":
        from llama_index.llms.openai import OpenAI as LIOpenAI
        _li_llm = LIOpenAI(model=model)
    else:
        from llama_index.llms.anthropic import Anthropic as LIAnthropic
        _li_llm = LIAnthropic(model=model)

    # Build A1 orchestrator
    print("Creating A1 orchestrator...")
    from biomni.agent.a1 import A1

    a1 = A1(
        path=str(agentic_root / ".biomni_data"),
        llm=model,
        source=source,
        use_tool_retriever=False,
        expected_data_lake_files=[],
    )

    # Register tools — clear default biomni tools to keep system prompt compact
    print("Registering tools...")
    a1.module2api = {}

    _register_tool(a1, _build_read_file_fn(_REPO_ROOT), _READ_FILE_SCHEMA)
    _register_tool(a1, _build_list_directory_fn(_REPO_ROOT), _LIST_DIRECTORY_SCHEMA)
    _register_tool(a1, _build_search_manuscript_fn(manuscript_index, _li_llm), _SEARCH_MANUSCRIPT_SCHEMA)

    a1.configure()
    a1._inject_custom_functions_to_repl()

    # Prepend framework overview to A1's system prompt
    prefix = _build_system_prompt_prefix(docs_dir)
    a1.system_prompt = prefix + a1.system_prompt

    print("A1 agent ready!")
    return a1

