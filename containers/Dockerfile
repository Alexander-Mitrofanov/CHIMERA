# syntax=docker/dockerfile:1

FROM python:3.11.15-slim-bookworm AS builder

ENV PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PYTHONDONTWRITEBYTECODE=1

WORKDIR /build

COPY pyproject.toml README.md LICENSE ./
COPY src/ ./src/
COPY schemas/ ./schemas/

RUN python -m pip install --no-cache-dir "build==1.2.2" \
    && python -m build --wheel --outdir /wheelhouse


FROM python:3.11.15-slim-bookworm AS runtime

ARG CHIMERA_VERSION=1.0.0
ARG VCS_REF=unknown

LABEL org.opencontainers.image.title="CHIMERA" \
      org.opencontainers.image.description="Leakage-aware virus-versus-host benchmark generator" \
      org.opencontainers.image.version="${CHIMERA_VERSION}" \
      org.opencontainers.image.revision="${VCS_REF}" \
      org.opencontainers.image.source="https://github.com/Alexander-Mitrofanov/CHIMERA" \
      org.opencontainers.image.documentation="https://github.com/Alexander-Mitrofanov/CHIMERA#readme" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.vendor="CHIMERA contributors"

ENV PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

RUN groupadd --gid 10001 chimera \
    && useradd --uid 10001 --gid chimera --create-home \
        --home-dir /home/chimera --shell /usr/sbin/nologin chimera \
    && install -d --owner chimera --group chimera /data

COPY --from=builder /wheelhouse/ /wheelhouse/

RUN python -m pip install --no-cache-dir --upgrade "pip>=26.1.2,<27" \
    && python -m pip install --no-cache-dir /wheelhouse/*.whl \
    && python -c "from chimera.schema_resources import validate_packaged_schemas; validate_packaged_schemas()" \
    && chimera --version \
    && rm -rf /wheelhouse

USER chimera
WORKDIR /data

HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD ["chimera", "--version"]

ENTRYPOINT ["chimera"]
CMD ["--help"]
