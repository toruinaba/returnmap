# gfortran + Python, with the project's dependencies (returnmap included) in
# /opt/venv.  The project itself is bind-mounted at /workspace at run time, so
# the image is rebuilt only when pyproject.toml changes.
#
#   make docker-image
#   make docker-check

FROM python:3.12-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        gfortran git make liblapack-dev libblas-dev \
    && rm -rf /var/lib/apt/lists/*

# Outside /workspace: that path is bind-mounted, which would hide the venv.
ENV UV_PROJECT_ENVIRONMENT=/opt/venv
ENV PATH=/opt/venv/bin:$PATH

WORKDIR /build
COPY pyproject.toml ./
RUN pip install --no-cache-dir uv \
 && uv sync --no-install-project --extra dev --extra fortran

WORKDIR /workspace
