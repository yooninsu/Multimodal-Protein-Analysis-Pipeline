FROM python:3.9-slim

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

# 시스템 패키지 설치
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    build-essential \
    libffi-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Python 의존성 설치
COPY requirements.txt .
RUN pip install --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# ColabFold 설치
RUN pip install colabfold[alphafold]

# 파일들 복사
COPY analyze_structure.py .
COPY run_pipeline.sh .

RUN chmod +x analyze_structure.py run_pipeline.sh

# 디렉토리 생성
RUN mkdir -p /app/input /app/output

ENTRYPOINT ["/app/run_pipeline.sh"]
CMD ["full"]