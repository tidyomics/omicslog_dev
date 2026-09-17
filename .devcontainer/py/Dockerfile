# Pin Python version here
FROM python:3.11-slim

ENV DEBIAN_FRONTEND=noninteractive

# System deps for common scientific wheels / builds
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    curl \
    build-essential \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# Create a non-root user for VS Code
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=1000
RUN groupadd --gid $USER_GID $USERNAME \
  && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
  && mkdir -p /workspaces \
  && chown -R $USERNAME:$USERNAME /workspaces

# Poetry install (pin version for reproducibility)
ENV POETRY_VERSION=1.8.3
RUN pip install --no-cache-dir "poetry==${POETRY_VERSION}"

# Poetry settings: keep venv in project for deterministic VS Code interpreter selection
ENV POETRY_VIRTUALENVS_IN_PROJECT=true \
    POETRY_NO_INTERACTION=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

USER $USERNAME
WORKDIR /workspaces/project

