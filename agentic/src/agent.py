"""
GeneRNBI orchestrating agent — powered by Biomni A1.

Two tools:
  - search_docs:   RAG over docs/source/ (API spec, datasets, integration guide)
  - search_viash:  Live page fetch from viash.io/guide + /reference (no pre-indexing)
"""

import os
from pathlib import Path
from typing import Optional

from dotenv import load_dotenv

from config import LLM_PROVIDER, MODEL_ID
from rag_builder import (
    _missing_key_error,
    get_or_build_docs_index,
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


_SEARCH_VIASH_SCHEMA = {
    "name": "search_viash",
    "description": (
        "Search the official Viash documentation (viash.io/guide and viash.io/reference) "
        "for information about Viash component configuration, Docker engine setup, "
        "Nextflow runner, argument types, resources, test benches, and CLI commands. "
        "Use this tool when the user has a Viash, Docker, or Nextflow related question "
        "during integration or troubleshooting — especially for config.vsh.yaml syntax, "
        "engine/runner setup, or Viash CLI usage."
    ),
    "required_parameters": [
        {
            "name": "query",
            "type": "str",
            "description": "Question about Viash, Docker, or Nextflow.",
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

def _build_docs_search_fn(index):
    """Direct retriever for docs — returns raw chunks without LLM synthesis."""
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


def _build_viash_live_search_fn():
    """Live search on viash.io — fetches pages on the fly, no pre-indexing.

    On each call:
      1. Fetches the viash.io sitemap (cached in memory after first call).
      2. Scores guide/reference URLs by keyword overlap with the query.
      3. Fetches and returns the text of the top matching pages.
    """
    import re
    import requests
    from bs4 import BeautifulSoup

    _cache = {"urls": None}

    def _get_urls():
        if _cache["urls"] is not None:
            return _cache["urls"]
        try:
            r = requests.get("https://viash.io/sitemap.xml", timeout=15)
            r.raise_for_status()
            all_urls = re.findall(r"<loc>(.*?)</loc>", r.text)
            urls = [
                u for u in all_urls
                if any(x in u for x in ["/guide/", "/reference/", "/faq"])
                and "/deprecated/" not in u
            ]
            _cache["urls"] = urls
            return urls
        except Exception as e:
            return []

    def _fetch_page(session, url):
        try:
            resp = session.get(url, timeout=10)
            if not resp.ok:
                return None
            soup = BeautifulSoup(resp.text, "html.parser")
            for tag in soup.find_all(["nav", "header", "footer", "aside", "script", "style"]):
                tag.decompose()
            main = (
                soup.find("main")
                or soup.find("article")
                or soup.find("div", class_=re.compile(r"content|prose|doc"))
            )
            text = (main or soup).get_text(separator="\n", strip=True)
            return text if len(text) >= 100 else None
        except Exception:
            return None

    def _search(query: str) -> str:
        try:
            keywords = [w for w in re.split(r"\W+", query.lower()) if len(w) > 2]
            urls = _get_urls()
            if not urls:
                return "Could not fetch viash.io sitemap. Check your internet connection."

            scored = sorted(
                urls,
                key=lambda u: sum(1 for kw in keywords if kw in u.lower()),
                reverse=True,
            )

            session = requests.Session()
            session.headers["User-Agent"] = "Mozilla/5.0 (viash-docs-searcher)"
            parts = []
            for url in scored:
                text = _fetch_page(session, url)
                if text:
                    parts.append(f"[Source: {url}]\n{text[:3000]}")
                if len(parts) >= 3:
                    break

            result = "\n\n---\n\n".join(parts) if parts else "No relevant viash.io pages found."
            print(result)
            return result
        except Exception as e:
            err = f"Error searching viash.io: {e}"
            print(err)
            return err

    return _search



# ─── System prompt ────────────────────────────────────────────────────────────

_SYSTEM_PROMPT_PREFIX = """\
You are the GeneRNBI Assistant — expert on the GeneRNBI GRN benchmarking framework.

══════════════════════════════════════════════════════
ABSOLUTE RULES — NEVER BREAK THESE
══════════════════════════════════════════════════════

RULE 1 — TOOLS FIRST, ALWAYS.
  Before writing ANY code, file, or answer, you MUST call search_docs.
  If the question involves Viash, Docker, or Nextflow, also call search_viash.
  Only use what the tools return. Never answer from memory.

RULE 2 — IMPLEMENTATION TASKS: CREATE REAL FILES ON DISK.
  When asked to write a method, you MUST:
    a) Call search_docs first to get the exact file conventions.
    b) Write script.py to disk using Python open() or pathlib.
    c) Write config.vsh.yaml to disk using Python open() or pathlib.
    d) Run viash commands using subprocess.run() — NEVER os.system().
       ALWAYS capture output so you can see errors:
         import subprocess
         result = subprocess.run(
             ["viash", "test", "path/to/config.vsh.yaml"],
             capture_output=True, text=True
         )
         print(result.stdout)
         print(result.stderr)
       Read the output. If it contains "ERROR" or "failed", report the error.
       NEVER claim success if the output contains an error message.
  Do NOT just print code. Do NOT just describe steps. CREATE the files.

RULE 3 — USE EXACTLY THE GENERNIBI FORMAT (from search_docs output).
  The correct config.vsh.yaml format uses:
    __merge__, name, namespace, info, resources, engines, runners
  NOT "version/components" — that is wrong and will fail.

RULE 4 — NEVER SKIP STEPS.
  If you are asked to "write and run" something, you must do both:
  create the files AND execute viash test.

RULE 5 — REPORT FAILURES HONESTLY.
  If a command fails (Docker not running, missing file, import error, etc.),
  say exactly what failed and why. Never claim "tested successfully" if there
  was an error in the output.

══════════════════════════════════════════════════════
KNOWN FILE TEMPLATES (use these as starting points)
══════════════════════════════════════════════════════

## config.vsh.yaml (GeneRNBI method format):
```yaml
__merge__: /src/api/comp_method.yaml

name: method_name
namespace: "grn_methods"
info:
  label: "Method Name"
  summary: "Short description"
resources:
  - type: python_script
    path: script.py
engines:
  - type: docker
    image: ghcr.io/openproblems-bio/base_python:1.0.4
    __merge__: /src/api/base_requirements.yaml
    setup:
      - type: python
        packages: [ numpy, pandas, scipy, anndata ]
runners:
  - type: executable
  - type: nextflow
    directives:
      label: [midtime, midmem, midcpu]
```

## script.py (GeneRNBI method format):
```python
import sys
import anndata as ad
import numpy as np
import pandas as pd

## VIASH START
par = {
    'rna': 'resources_test/grn_benchmark/rna_op.h5ad',
    'tf_all': 'resources_test/grn_benchmark/prior/tf_all.csv',
    'prediction': 'output/prediction.h5ad'
}
## VIASH END

rna = ad.read_h5ad(par['rna'])
tf_all = np.loadtxt(par['tf_all'], dtype=str)

# <<< YOUR METHOD HERE >>>

net['weight'] = net['weight'].astype(str)
output = ad.AnnData(
    X=None,
    uns={
        'method_id': 'method_name',
        'dataset_id': rna.uns.get('dataset_id', 'unknown'),
        'prediction': net[['source', 'target', 'weight']]
    }
)
output.write(par['prediction'])
```

## viash test command:
```bash
viash test src/methods/your_method/config.vsh.yaml
```

══════════════════════════════════════════════════════
TOOL USAGE
══════════════════════════════════════════════════════

search_docs  — ALWAYS call this first for any GeneRNBI question:
  • File structure, config.vsh.yaml format, script.py conventions
  • Datasets, metrics, API formats
  Usage:  result = search_docs('your question'); print(result)

search_viash  — Call this for Viash/Docker/Nextflow details:
  • config.vsh.yaml engine/runner syntax
  • Docker setup, apt/python packages
  • viash CLI commands (viash test, viash run, viash build)
  Usage:  result = search_viash('your question'); print(result)

"""


# ─── Main initializer ─────────────────────────────────────────────────────────

def initialize_agent(
    docs_dir: Optional[Path] = None,
    use_cache: bool = True,
    model_id: str = None,
    llm_provider: str = None,
):
    """Build docs index and return an A1 orchestrator with two tools."""
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
    repo_root = agentic_root.parent
    if docs_dir is None:
        docs_dir = repo_root / "docs" / "source"

    print("Building docs index...")
    docs_index = get_or_build_docs_index(docs_dir, use_cache=use_cache)

    # LlamaIndex LLM no longer needed — all tools use direct retrieval, no synthesis
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

    search_docs_fn = _build_docs_search_fn(docs_index)
    search_docs_fn.__doc__ = (
        "Search GeneRNBI documentation for integration guides, API formats, "
        "dataset info, and framework usage.\n\n"
        "Args:\n    query (str): Any technical or integration question.\n\n"
        "Returns:\n    str: Answer from the documentation."
    )
    search_docs_fn.__name__ = "search_docs"

    search_viash_fn = _build_viash_live_search_fn()
    search_viash_fn.__doc__ = (
        "Search viash.io documentation for Viash config syntax, Docker engine setup, "
        "Nextflow runner configuration, argument types, and Viash CLI usage.\n\n"
        "Args:\n    query (str): Question about Viash, Docker, or Nextflow.\n\n"
        "Returns:\n    str: Answer from viash.io docs."
    )
    search_viash_fn.__name__ = "search_viash"

    search_manuscript_fn = None  # removed

    _register_tool(a1, search_docs_fn, _SEARCH_DOCS_SCHEMA)
    _register_tool(a1, search_viash_fn, _SEARCH_VIASH_SCHEMA)

    a1.configure()
    a1._inject_custom_functions_to_repl()

    a1.system_prompt = _SYSTEM_PROMPT_PREFIX + a1.system_prompt

    print("A1 agent ready!")
    return a1


