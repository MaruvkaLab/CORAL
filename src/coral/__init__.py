"""CORAL: Multi-Species Mutation Extraction Pipeline"""

__all__ = [
    "MutationExtractionPipeline",
    "MultiSpeciesMutationPipeline",
]


def __getattr__(name):
    """Import the pipelines on first use."""
    if name in __all__:
        from . import pipeline
        return getattr(pipeline, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

