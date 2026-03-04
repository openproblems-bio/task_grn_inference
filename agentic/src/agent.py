"""
GeneRNBI orchestrating agent — powered by Biomni A1.

Two tools:
  - search_docs:   RAG over docs/source/ (API spec, datasets, integration guide)
  - search_viash:  Live page fetch from viash.io/guide + /reference (no pre-indexing)

Know-how documents in agentic/know_how/ are loaded into the agent's system prompt
via Biomni's KnowHowLoader, following the same pattern as the upstream Biomni project.
"""

import os
import re
import requests
from pathlib import Path
from typing import Optional

from dotenv import load_dotenv

from config import LLM_PROVIDER, MODEL_ID
from rag_builder import (
    _missing_key_error,
    get_or_build_docs_index,
)

_SOURCE_MAP = {"openai": "OpenAI", "anthropic": "Anthropic"}


# ─── Tool implementations ─────────────────────────────────────────────────────

def _build_docs_search_fn(index):
    """Return a search_docs function over the LlamaIndex docs index."""
    retriever = index.as_retriever(similarity_top_k=10)

    def search_docs(query: str) -> str:
        """Search the GeneRNBI documentation for technical information about the framework.

        Use this tool FIRST for any question about:
        - How to integrate a new GRN inference method or metric
        - config.vsh.yaml and script.py file format and conventions
        - Available datasets, their format, and how to download them
        - Output format for inferred GRN networks (AnnData, source/target/weight columns)
        - Running evaluation metrics and interpreting scores
        - Any GeneRNBI API or pipeline question

        Args:
            query (str): Natural language question about GeneRNBI usage, integration, or API.

        Returns:
            str: Relevant documentation excerpts.
        """
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

    return search_docs


def _build_viash_live_search_fn():
    """Return a search_viash function that fetches viash.io pages on demand."""
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
        except Exception:
            return []

    def _fetch_page(session, url):
        from bs4 import BeautifulSoup
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

    def search_viash(query: str) -> str:
        """Search the official Viash documentation (viash.io) for Viash/Docker/Nextflow help.

        Use this tool when the user has a question about:
        - config.vsh.yaml syntax: engines, runners, resources, arguments
        - Docker engine setup: image, setup, aptRequirements, pythonRequirements
        - Nextflow runner configuration and directives
        - Viash CLI commands: viash test, viash run, viash build, viash ns build
        - Troubleshooting Viash, Docker image builds, or Nextflow pipeline errors

        This tool fetches pages live from viash.io — no pre-built index needed.

        Args:
            query (str): Question about Viash, Docker, or Nextflow configuration.

        Returns:
            str: Relevant content from viash.io documentation pages.
        """
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

    return search_viash


# ─── Main initializer ─────────────────────────────────────────────────────────

def initialize_agent(
    docs_dir: Optional[Path] = None,
    use_cache: bool = True,
    model_id: str = None,
    llm_provider: str = None,
):
    """Build docs index and return an A1 orchestrator with two tools and know-how docs."""
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

    print("Creating A1 orchestrator...")
    from biomni.agent.a1 import A1
    from biomni.know_how import KnowHowLoader

    # Load domain-specific know-how documents from agentic/know_how/
    know_how_dir = agentic_root / "know_how"
    know_how_loader = KnowHowLoader(know_how_dir=str(know_how_dir))
    print(f"Loaded {len(know_how_loader.documents)} know-how documents from {know_how_dir}")

    a1 = A1(
        path=str(agentic_root / ".biomni_data"),
        llm=model,
        source=source,
        use_tool_retriever=False,
        expected_data_lake_files=[],
    )

    # Replace default know-how loader with our domain-specific one
    # Must be done before configure() is called (add_tool triggers it)
    a1.know_how_loader = know_how_loader

    # Register tools using add_tool() — auto-generates schema from docstring
    # and triggers configure() (which injects know-how into system prompt)
    print("Registering tools...")
    a1.add_tool(_build_docs_search_fn(docs_index))
    a1.add_tool(_build_viash_live_search_fn())

    print("A1 agent ready!")
    return a1
