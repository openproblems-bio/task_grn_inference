"""
Build and manage RAG indices.

  - docs index:       indexes agentic/docs/ (markdown documentation)
  - manuscript index: indexes agentic/data/ (PDF manuscript)
"""

import sys
from pathlib import Path


def _missing_key_error():
    print(
        "\n"
        "╔══════════════════════════════════════════════════════════════╗\n"
        "║  ✗  API key not found                                        ║\n"
        "╠══════════════════════════════════════════════════════════════╣\n"
        "║  Set your key using ONE of these methods:                    ║\n"
        "║                                                              ║\n"
        "║  1. Inline (one-off run):                                    ║\n"
        "║     OPENAI_API_KEY=sk-... bash run.sh                        ║\n"
        "║                                                              ║\n"
        "║  2. Export for the session:                                  ║\n"
        "║     export OPENAI_API_KEY=sk-...                             ║\n"
        "║                                                              ║\n"
        "║  3. Permanent (add to ~/.bashrc or ~/.zshrc):                ║\n"
        "║     echo 'export OPENAI_API_KEY=sk-...' >> ~/.bashrc         ║\n"
        "║                                                              ║\n"
        "║  4. .env file (copy .env.template → .env, fill in key)       ║\n"
        "║                                                              ║\n"
        "║  Get your API key:                                           ║\n"
        "║  • OpenAI:    https://platform.openai.com/api-keys           ║\n"
        "║  • Anthropic: https://console.anthropic.com/settings/keys    ║\n"
        "╚══════════════════════════════════════════════════════════════╝\n",
        file=sys.stderr,
    )
    sys.exit(1)


from llama_index.core import (
    VectorStoreIndex,
    SimpleDirectoryReader,
    StorageContext,
    load_index_from_storage,
)


def _build_index(source_dir: Path, persist_dir: Path, label: str, **reader_kwargs) -> VectorStoreIndex:
    """Shared helper: load from cache or build from source_dir."""
    if not source_dir.exists():
        raise ValueError(f"{label} folder not found at {source_dir}")

    if persist_dir and persist_dir.exists():
        print(f"Loading cached {label} index from {persist_dir}")
        storage_context = StorageContext.from_defaults(persist_dir=str(persist_dir))
        return load_index_from_storage(storage_context)

    print(f"Building {label} index from {source_dir}")
    documents = SimpleDirectoryReader(str(source_dir), **reader_kwargs).load_data()
    print(f"  Loaded {len(documents)} document(s)")

    try:
        index = VectorStoreIndex.from_documents(documents, show_progress=True)
    except ValueError as e:
        if "api key" in str(e).lower() or "openai" in str(e).lower():
            _missing_key_error()
        raise

    if persist_dir:
        persist_dir.mkdir(parents=True, exist_ok=True)
        index.storage_context.persist(persist_dir=str(persist_dir))
        print(f"  {label.capitalize()} index persisted to {persist_dir}")

    return index


def get_or_build_docs_index(docs_dir: Path, use_cache: bool = True) -> VectorStoreIndex:
    """Build or load the docs index (agentic/docs/)."""
    persist_dir = Path(__file__).parent.parent / ".index_cache" / "docs" if use_cache else None
    return _build_index(docs_dir, persist_dir, label="docs", recursive=True)


def get_or_build_manuscript_index(data_dir: Path, use_cache: bool = True) -> VectorStoreIndex:
    """Build or load the manuscript index (agentic/data/ — PDF)."""
    persist_dir = Path(__file__).parent.parent / ".index_cache" / "manuscript" if use_cache else None
    return _build_index(data_dir, persist_dir, label="manuscript", recursive=False)
