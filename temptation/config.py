"""Thresholds and tunables in one place.

No magic numbers should appear anywhere else in the package. Defaults below
match the current CLI defaults in measure_nerve.py exactly, so building a
SegmentationConfig with no overrides reproduces today's behavior.
"""

import dataclasses
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass(frozen=True)
class SegmentationConfig:
    myelin_val: int = 64
    axoplasm_val: int = 192
    mito_val: int = 128

    mode: str = "auto"                          # "normal" | "pathological" | "auto"
    myelin_threshold_px: int = 200

    smoothing_radius_px: int = 1
    min_axon_area_px: int = 200
    min_myelin_area_px: int = 300

    watershed_mode: str = "weighted"             # "simple" | "weighted"
    watershed_weight: str = "radius"              # "radius" | "area"
    watershed_compactness: float = 0.001
    watershed_beta: float = 1.0

    assign_detached_myelin: str = "none"          # "none" | "nearest"

    # Phase 2b (gated): "legacy" reproduces the pre-fix mitochondrial-rim
    # bug (F1) exactly; "fill" is not implemented until Phase 2b lands.
    mito_hole_handling: str = "legacy"

    # Phase 3: mitochondria-to-axon assignment strategy.
    mito_assignment: str = "legacy"


@dataclass(frozen=True)
class QCThresholds:
    g_ratio_min: float = 0.3
    g_ratio_max: float = 0.95
    axon_circularity_min: float = 0.4
    min_myelin_area_um2: Optional[float] = None
    min_axon_area_um2: Optional[float] = None
    min_fiber_area_um2: Optional[float] = None
    min_mito_count_for_clustering: int = 3


@dataclass(frozen=True)
class ComponentSpec:
    metric: str
    op: str                         # ">" or "<"
    threshold: Optional[float]
    weight: float
    enabled: bool = True


@dataclass(frozen=True)
class PathologyScoreConfig:
    enabled: bool = False
    version: str = "v0"
    components: dict = field(default_factory=dict)


@dataclass(frozen=True)
class AnalysisConfig:
    segmentation: SegmentationConfig = field(default_factory=SegmentationConfig)
    qc: QCThresholds = field(default_factory=QCThresholds)
    pathology_score: PathologyScoreConfig = field(default_factory=PathologyScoreConfig)
    pixel_size_um: Optional[float] = None
    group: Optional[str] = None

    def to_dict(self) -> dict:
        return dataclasses.asdict(self)


def load_config(path) -> AnalysisConfig:
    """Load an AnalysisConfig from a JSON file (YAML support added when
    pyyaml is available; not a dependency of this package)."""
    path = Path(path)
    text = path.read_text()
    if path.suffix in (".yaml", ".yml"):
        try:
            import yaml
        except ImportError as exc:
            raise ImportError(
                "pyyaml is required to load .yaml config files; "
                "pip install pyyaml, or use a .json config instead."
            ) from exc
        data = yaml.safe_load(text)
    else:
        data = json.loads(text)

    seg = SegmentationConfig(**data.get("segmentation", {}))
    qc = QCThresholds(**data.get("qc_thresholds", data.get("qc", {})))
    score_data = data.get("pathology_score", {})
    score = PathologyScoreConfig(
        enabled=score_data.get("enabled", False),
        version=score_data.get("version", "v0"),
        components=score_data.get("components", {}),
    )
    return AnalysisConfig(
        segmentation=seg,
        qc=qc,
        pathology_score=score,
        pixel_size_um=data.get("pixel_size_um"),
        group=data.get("group"),
    )
