"""
GeneRNBI orchestrating agent — powered by Biomni A1.

Two tools (registered via a1.add_tool() — standard Biomni API):
  - search_docs:   RAG over docs/source/ (API spec, datasets, integration guide)
  - search_viash:  Live page fetch from viash.io/guide + /reference (no pre-indexing)

Know-how documents in agentic/know_how/ are injected into A1's know_how_loader
so configure() includes them in the system prompt automatically.
"""

import os
import re
from pathlib import Path
from typing import Optional

from dotenv import load_dotenv

from config import LLM_PROVIDER, MODEL_ID
from rag_builder import (
    _missing_key_error,
    get_or_build_docs_index,
)

_SOURCE_MAP = {"openai": "OpenAI", "anthropic": "Anthropic"}

# Module-level state — set by initialize_agent() before add_tool() calls so
# that inspect.getsource() on the top-level functions works cleanly.
_DOCS_RETRIEVER = None
_VIASH_URL_CACHE: dict = {"urls": None}


# ─── Tool functions (module-level so inspect.getsource captures them correctly) ──

def search_docs(query: str) -> str:
    """Search the GeneRNBI documentation for information about integrating new
    methods, metrics, or datasets; the API format for Viash components
    (config.vsh.yaml, script.py structure, argument types, file formats);
    dataset descriptions; and any other framework usage or integration question.
    Use this tool FIRST for all technical and integration questions.

    Args:
        query (str): Any question about framework usage, integration, or API structure.

    Returns:
        str: Relevant documentation content.
    """
    global _DOCS_RETRIEVER
    try:
        nodes = _DOCS_RETRIEVER.retrieve(query)
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


def search_viash(query: str) -> str:
    """Search the official Viash documentation (viash.io/guide and viash.io/reference)
    for information about Viash component configuration, Docker engine setup,
    Nextflow runner, argument types, resources, test benches, and CLI commands.
    Use this tool when the user has a Viash, Docker, or Nextflow related question
    during integration or troubleshooting — especially for config.vsh.yaml syntax,
    engine/runner setup, or Viash CLI usage.

    Args:
        query (str): Question about Viash, Docker, or Nextflow.

    Returns:
        str: Relevant content fetched live from viash.io.
    """
    import requests
    from bs4 import BeautifulSoup

    def _get_urls():
        if _VIASH_URL_CACHE["urls"] is not None:
            return _VIASH_URL_CACHE["urls"]
        try:
            r = requests.get("https://viash.io/sitemap.xml", timeout=15)
            r.raise_for_status()
            all_urls = re.findall(r"<loc>(.*?)</loc>", r.text)
            urls = [
                u for u in all_urls
                if any(x in u for x in ["/guide/", "/reference/", "/faq"])
                and "/deprecated/" not in u
            ]
            _VIASH_URL_CACHE["urls"] = urls
            return urls
        except Exception:
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


# ─── Main initializer ─────────────────────────────────────────────────────────

def initialize_agent(
    docs_dir: Optional[Path] = None,
    use_cache: bool = True,
    model_id: str = None,
    llm_provider: str = None,
):
    """Build docs index and return an A1 orchestrator with two tools and know-how docs."""
    global _DOCS_RETRIEVER

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
    _DOCS_RETRIEVER = docs_index.as_retriever(similarity_top_k=10)

    print("Creating A1 orchestrator...")
    from biomni.agent.a1 import A1

    a1 = A1(
        path=str(agentic_root / ".biomni_data"),
        llm=model,
        source=source,
        use_tool_retriever=False,
        expected_data_lake_files=[],
    )

    # Inject know-how docs into A1's know_how_loader so configure() includes them
    know_how_dir = agentic_root / "know_how"
    if know_how_dir.exists():
        n_loaded = 0
        for md_path in sorted(know_how_dir.glob("*.md")):
            try:
                content = md_path.read_text(encoding="utf-8")
                desc = ""
                for line in content.splitlines():
                    if line.startswith("**Short Description**"):
                        desc = line.split(":", 1)[-1].strip().strip("*").strip()
                        break
                a1.know_how_loader.add_custom_document(
                    doc_id=f"grn_{md_path.stem}",
                    name=md_path.stem.replace("_", " ").title(),
                    description=desc or md_path.stem,
                    content=content,
                )
                n_loaded += 1
            except Exception as e:
                print(f"Warning: could not load {md_path.name}: {e}")
        print(f"Injected {n_loaded} know-how doc(s) into A1")

    # Register tools via the standard Biomni add_tool() API.
    # Module-level functions ensure inspect.getsource() captures them cleanly.
    print("Registering tools...")
    a1.add_tool(search_docs)
    a1.add_tool(search_viash)

    print("A1 agent ready!")
    return a1
