
## BUILD USING: docker build -t knitwork:latest .
## RUN USING: 

# docker run -v $(pwd):/app/data -e NEO4J_LOCATION=$NEO4J_LOCATION -e NEO4J_USERNAME=$NEO4J_USERNAME -e NEO4J_PASSWORD=$NEO4J_PASSWORD knitwork:latest 

FROM python:3.11-slim

WORKDIR /app

RUN apt-get update \
 && apt-get install -y --no-install-recommends libxrender1 libxext6 libsm6 libexpat1 \
 && rm -rf /var/lib/apt/lists/*

COPY pyproject.toml README.md .

RUN pip install --no-cache-dir .

COPY knitwork knitwork
COPY run_knitwork.sh .

RUN chmod +x run_knitwork.sh

CMD ["bash", "run_knitwork.sh"]
