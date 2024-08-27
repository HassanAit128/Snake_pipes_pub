FROM python:3.9-slim

# Update and install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    libssl-dev \
    libffi-dev \
    python3-dev \
    python3-pip \
    && apt-get clean

# Install MACS2
RUN pip install MACS2

# Set environment variables
ENV LC_ALL=C
ENV LANG=C