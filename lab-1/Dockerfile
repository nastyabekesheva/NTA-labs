FROM python:3.11-slim

WORKDIR /src
COPY src/. .

RUN apt update && apt install -y gcc
RUN apt-get update

RUN pip install --upgrade pip
RUN pip install numpy
RUN pip install sympy

CMD ["python", "script.py"]
