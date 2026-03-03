"""
Integration sub-agent for the GRN benchmark.

Responsibilities:
  - Answer questions about how to integrate new methods, metrics, or datasets
  - Read and explain existing source code in src/metrics/ and src/methods/
  - Fetch and interpret Viash documentation (viash.io)
  - Troubleshoot pipeline and configuration errors

Tools available:
  search_docs        — RAG over agentic/docs/
  list_source_files  — list files under src/metrics/ or src/methods/
  read_source_file   — read a specific source file
  fetch_viash_docs   — fetch a viash.io documentation page
"""

import re
from html.parser import HTMLParser
from pathlib import Path
from urllib.request import urlopen
from urllib.error import URLError

import requests
from llama_index.core import VectorStoreIndex
from llama_index.core.agent.workflow import ReActAgent
from llama_index.core.tools import FunctionTool

from config import LLM_PROVIDER, MODEL_ID
from rag_builder import get_or_build_docs_index

# Root of the task_grn_inference source tree (two levels above agentic/)
_SRC_ROOT = Path(__file__).parent.parent.parent / "src"
_ALLOWED_DIRS = (_SRC_ROOT / "metrics", _SRC_ROOT / "methods")


# ---------------------------------------------------------------------------
# Tool implementations
# ---------------------------------------------------------------------------

class _HTMLTextExtractor(HTMLParser):
    def __init__(self):
        super().__init__()
        self._skip_tags = {"script", "style", "nav", "footer", "header"}
        self._current_skip = 0
        self.text_parts = []

    def handle_starttag(self, tag, attrs):
        if tag in self._skip_tags:
            self._current_skip += 1

    def handle_endtag(self, tag):
        if tag in self._skip_tags and self._current_skip:
            self._current_skip -= 1

    def handle_data(self, data):
        if not self._current_skip:
            stripped = data.strip()
            if stripped:
                self.text_parts.append(stripped)


def _fetch_url_text(url: str, max_chars: int = 8000) -> str:
    try:
        resp = requests.get(url, timeout=10, headers={"User-Agent": "Mozilla/5.0"})
        resp.raise_for_status()
        parser = _HTMLTextExtractor()
        parser.feed(resp.text)
        text = "\n".join(parser.text_parts)
        # collapse excessive blank lines
        text = re.sub(r"\n{3,}", "\n\n", text)
        return text[:max_chars]
    except Exception as e:
        return f"Error fetching {url}: {e}"


def _resolve_path(filepath: str) -> Path | None:
    """Resolve filepath to an absolute path, restricted to allowed directories."""
    p = Path(filepath)
    if not p.is_absolute():
        # try relative to each allowed dir
        for allowed in _ALLOWED_DIRS:
            candidate = (allowed / p).resolve()
            if candidate.exists():
                return candidate
        # try relative to _SRC_ROOT
        candidate = (_SRC_ROOT / p).resolve()
        if candidate.exists():
            return candidate
        return None
    p = p.resolve()
    if any(str(p).startswith(str(d)) for d in _ALLOWED_DIRS):
        return p
    return None


def list_source_files(directory: str = "") -> str:
    """List source files available under src/metrics/ or src/methods/.

    Args:
        directory: Subdirectory to list, e.g. 'metrics', 'methods', 'metrics/regression'.
                   Leave empty to list both top-level directories.
    """
    if directory:
        target = (_SRC_ROOT / directory).resolve()
        if not target.exists():
            return f"Directory not found: {directory}"
        dirs_to_list = [target]
    else:
        dirs_to_list = list(_ALLOWED_DIRS)

    lines = []
    for d in dirs_to_list:
        lines.append(f"\n📁 {d.relative_to(_SRC_ROOT.parent)}/")
        for f in sorted(d.rglob("*")):
            if f.is_file() and f.suffix in {".py", ".yaml", ".yml", ".md", ".txt"}:
                lines.append(f"   {f.relative_to(_SRC_ROOT)}")
    return "\n".join(lines) if lines else "No files found."


def read_source_file(filepath: str) -> str:
    """Read the contents of a source file from src/metrics/ or src/methods/.

    Args:
        filepath: Path relative to src/, e.g. 'metrics/regression/script.py'
                  or 'methods/grnboost/config.vsh.yaml'.
    """
    resolved = _resolve_path(filepath)
    if resolved is None:
        return (
            f"File not found or not accessible: '{filepath}'.\n"
            f"Use list_source_files() to see available files. "
            f"Only files under src/metrics/ and src/methods/ are accessible."
        )
    try:
        content = resolved.read_text(encoding="utf-8", errors="replace")
        rel = resolved.relative_to(_SRC_ROOT.parent)
        return f"# {rel}\n\n{content}"
    except Exception as e:
        return f"Error reading file: {e}"


def fetch_viash_docs(url: str = "https://viash.io/guide/") -> str:
    """Fetch a page from the Viash documentation site (viash.io).

    Use this to look up Viash component structure, config.vsh.yaml syntax,
    argument types, runner configuration, or packaging instructions.

    Args:
        url: A viash.io URL, e.g. 'https://viash.io/guide/' or
             'https://viash.io/reference/config/'.
    """
    if "viash.io" not in url:
        return "This tool only fetches pages from viash.io."
    return _fetch_url_text(url)


# ---------------------------------------------------------------------------
# Agent factory
# ---------------------------------------------------------------------------

def create_integration_agent(
    docs_index: VectorStoreIndex,
    model_id: str = None,
    llm_provider: str = None,
    verbose: bool = False,
) -> ReActAgent:
    """Create the integration sub-agent with all its tools."""
    model_id = model_id or MODEL_ID
    llm_provider = llm_provider or LLM_PROVIDER

    if llm_provider.lower() == "openai":
        from llama_index.llms.openai import OpenAI
        llm = OpenAI(model=model_id)
    elif llm_provider.lower() == "anthropic":
        from llama_index.llms.anthropic import Anthropic
        llm = Anthropic(model=model_id)
    else:
        raise ValueError(f"Unsupported LLM provider: {llm_provider}")

    docs_engine = docs_index.as_query_engine(similarity_top_k=5, llm=llm)

    def search_docs(query: str) -> str:
        """Search GeneRNBI/benchmark documentation for pipeline usage,
        dataset info, metric definitions, and integration guides.

        Args:
            query: Question or keywords about the benchmark framework.
        """
        try:
            return str(docs_engine.query(query))
        except Exception as e:
            return f"Error searching docs: {e}"

    tools = [
        FunctionTool.from_defaults(fn=search_docs, name="search_docs"),
        FunctionTool.from_defaults(fn=list_source_files, name="list_source_files"),
        FunctionTool.from_defaults(fn=read_source_file, name="read_source_file"),
        FunctionTool.from_defaults(fn=fetch_viash_docs, name="fetch_viash_docs"),
    ]

    system_prompt = (Path(__file__).parent / "instructions_integration.txt").read_text()
    return ReActAgent(llm=llm, tools=tools, system_prompt=system_prompt, verbose=verbose)


def initialize_integration_agent(
    docs_dir: Path = None,
    use_cache: bool = True,
    model_id: str = None,
    llm_provider: str = None,
) -> ReActAgent:
    """Build docs index and create the integration agent."""
    if docs_dir is None:
        docs_dir = Path(__file__).parent.parent / "docs"
    docs_index = get_or_build_docs_index(docs_dir, use_cache=use_cache)
    return create_integration_agent(docs_index, model_id=model_id, llm_provider=llm_provider)
