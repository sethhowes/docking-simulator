# Use the official Ubuntu 20.04 base image
FROM rayproject/ray:2.9.3-py39-cu118

# Set environment variable to non-interactive (this prevents some prompts)
ENV DEBIAN_FRONTEND=noninteractive
USER root

# Update package list and install required dependencies
RUN apt-get update \
    && apt-get install -y libboost-system-dev libboost-thread-dev \
    libboost-serialization-dev libboost-filesystem-dev \
    libboost-program-options-dev libboost-timer-dev \
    wget openbabel \
    git python3-pip

# Install Uni-Dock binary
RUN cd /opt && \
    wget https://github.com/dptech-corp/Uni-Dock/releases/download/1.0.0/unidock && \
    chmod +x unidock

# Ensure Uni-Dock binary is in path
ENV PATH="/opt:${PATH}"

# Install UniDockTools
RUN pip install git+https://github.com/dptech-corp/Uni-Dock.git#subdirectory=unidock_tools

# Install python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

WORKDIR /app/src/

# Copy source folder into app folder
COPY ./src /app/src

# Create an output directory
RUN mkdir ./data

# Run ray serve
CMD ["serve", "run", "-p", "1456", "docking-simulator.unidock:docking_app"]
