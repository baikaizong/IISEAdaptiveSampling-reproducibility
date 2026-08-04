"""Compress raster payloads inside EPS files while preserving vector content."""

from __future__ import annotations

import argparse
import base64
from io import BytesIO
from pathlib import Path
import re
import zlib

from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
FINAL = ROOT / "FinalVision"
TARGET_BYTES = 5_000_000
SAFETY_TARGET_BYTES = 4_800_000

IMAGE_BLOCK = re.compile(
    rb"/DataString\s+(?P<chunk>\d+)\s+string\s+def\s+"
    rb"(?P<width>\d+)\s+(?P<height>\d+)\s+8\s+"
    rb"(?P<matrix>\[[^\]]+\])\s+"
    rb"\{\s*currentfile\s+DataString\s+readhexstring\s+pop\s*\}\s+"
    rb"bind\s+false\s+3\s+colorimage\s*"
    rb"(?P<payload>[0-9A-Fa-f\s]+?)"
    rb"(?=\s*grestore)",
    re.DOTALL,
)


def ascii85(payload: bytes) -> bytes:
    """Encode a binary PostScript stream using wrapped Adobe ASCII85 syntax."""
    return base64.a85encode(payload, adobe=False, wrapcol=120) + b"~>\n"


def encoded_block(match: re.Match[bytes], mode: str, quality: int = 98) -> bytes:
    """Replace one Matplotlib hex RGB block with a filtered data source."""
    width = int(match.group("width"))
    height = int(match.group("height"))
    matrix = match.group("matrix")
    raw = bytes.fromhex(match.group("payload").decode("ascii"))
    expected = width * height * 3
    if len(raw) != expected:
        raise ValueError(
            f"Unexpected RGB payload: {len(raw)} bytes for {width}x{height}; expected {expected}"
        )

    if mode == "flate":
        compressed = zlib.compress(raw, level=9)
        decoder = b"/FlateDecode filter"
    elif mode == "jpeg":
        stream = BytesIO()
        Image.frombytes("RGB", (width, height), raw).save(
            stream,
            format="JPEG",
            quality=quality,
            subsampling=0,
            optimize=True,
            progressive=False,
        )
        compressed = stream.getvalue()
        decoder = b"/DCTDecode filter"
    else:
        raise ValueError(f"Unknown compression mode: {mode}")

    return b"\n".join(
        [
            f"/DataString {width * 3} string def".encode("ascii"),
            b"/DataSource currentfile /ASCII85Decode filter " + decoder + b" def",
            f"{width} {height} 8 ".encode("ascii") + matrix,
            b"{ DataSource DataString readstring pop } bind false 3 colorimage",
            ascii85(compressed),
        ]
    )


def rewrite(source: bytes, mode: str, quality: int = 98) -> tuple[bytes, int]:
    """Rewrite all recognized RGB blocks and return output plus block count."""
    count = 0

    def replace(match: re.Match[bytes]) -> bytes:
        nonlocal count
        count += 1
        return encoded_block(match, mode=mode, quality=quality)

    return IMAGE_BLOCK.sub(replace, source), count


def optimize(path: Path, output: Path) -> str:
    """Compress an EPS in place or to ``output`` without changing vector marks."""
    source = path.read_bytes()
    if len(source) <= TARGET_BYTES:
        if output != path:
            output.write_bytes(source)
        return f"keep {path.name}: {len(source) / 1_000_000:.2f} MB"

    candidate, count = rewrite(source, mode="flate")
    if count == 0:
        raise ValueError(f"No Matplotlib RGB image blocks found in {path}")

    method = "lossless Flate"
    if len(candidate) > SAFETY_TARGET_BYTES:
        for quality in (98, 97, 96, 95, 94, 92, 90):
            candidate, _ = rewrite(source, mode="jpeg", quality=quality)
            method = f"JPEG q={quality}, 4:4:4"
            if len(candidate) <= SAFETY_TARGET_BYTES:
                break

    if len(candidate) > TARGET_BYTES:
        raise ValueError(
            f"Cannot compress {path.name} below 5 MB without stronger quality reduction"
        )
    if not candidate.startswith(b"%!PS-Adobe"):
        raise ValueError(f"Invalid optimized EPS header for {path.name}")

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_bytes(candidate)
    temporary.replace(output)
    return (
        f"optimize {path.name}: {len(source) / 1_000_000:.2f} -> "
        f"{len(candidate) / 1_000_000:.2f} MB ({method}, {count} raster layer(s))"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*", type=Path)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()

    paths = args.paths or sorted(FINAL.glob("*.eps"))
    for path in paths:
        output = args.output_dir / path.name if args.output_dir else path
        print(optimize(path, output))


if __name__ == "__main__":
    main()
