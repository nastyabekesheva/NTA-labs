FROM python:3.11-slim

WORKDIR /

RUN apt update && apt install -y gcc
RUN apt-get update

RUN pip install --upgrade pip
RUN pip install numpy
RUN pip install sympy
RUN pip install math
RUN pip install sys

COPY script.py /
ENTRYPOINT ["python", "script.py"]
