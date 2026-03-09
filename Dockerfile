
# Official images are build using the GitHib CI 'latest' workflow,
# which runs automatically on commits to the main branch
# creating a new 'xchem/knitwork:latest' image.
#
# Otherwise, for local development you can build using:
#   docker build -t knitwork:latest .
# And ru using:
#   docker run -v $(pwd):/app/data -e NEO4J_LOCATION=$NEO4J_LOCATION -e NEO4J_USERNAME=user -e NEO4J_PASSWORD=pass knitwork:latest

FROM python:3.11-slim

WORKDIR /app

RUN chmod -R a+rw /app \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
      libxrender1 \
      libxext6 \
      libsm6 \
      libexpat1 \
    && rm -rf /var/lib/apt/lists/*

COPY pyproject.toml README.md ./

RUN pip install --no-cache-dir .

COPY knitwork knitwork
COPY run_knitwork.sh .

RUN chmod +x run_knitwork.sh

# Disable some advanced console (rich) features: -
ENV TERM=dumb
ENV TTY_COMPATIBLE=0
ENV TTY_INTERACTIVE=0

CMD ["bash", "run_knitwork.sh"]
