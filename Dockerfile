FROM python:3.11-slim

WORKDIR /app

# System dependencies for BioPython/pydna compilation
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
RUN pip install --no-cache-dir \
    biopython==1.84 \
    pydna==5.2.0 \
    primers==0.5.10 \
    pandas==2.2.2 \
    openpyxl==3.1.5 \
    numpy==2.0.0 \
    python-dotenv==1.0.1 \
    fastapi>=0.100.0 \
    "uvicorn[standard]>=0.23.0" \
    httpx \
    requests \
    pyparsing \
    pytest

COPY . .

# assets/ directory must contain pUC19.gb
# Place addgene-plasmid-50005-sequence-222046.gb as assets/pUC19.gb before building
RUN mkdir -p assets

EXPOSE 8001

CMD ["uvicorn", "api:app", "--host", "0.0.0.0", "--port", "8001"]
