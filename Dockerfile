FROM rayproject/ray-ml:2.6.3-gpu

RUN apt-get update && apt install -y cmake \
    libboost-system-dev libboost-thread-dev libboost-serialization-dev \
    libboost-filesystem-dev libboost-program-options-dev libboost-timer-dev \
    libsm6 \
    wget \
    openbabel \
    apt git python3-pip && \
    ln -s /usr/bin/python3 /usr/bin/python

# Install uni-dock
RUN cd opt && \
    wget https://github.com/dptech-corp/Uni-Dock/releases/download/1.0.0/unidock && \
    chmod +x unidock

# Ensure binaries are in path
ENV PATH="/opt:${PATH}"

# Install UniDockTools
RUN pip install git+https://github.com/dptech-corp/Uni-Dock.git#subdirectory=unidock_tools

# Install python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# # Copy source folder into app folder
# COPY ./src /app/src

# # Create an output directory
# RUN mkdir ./data

# # Move into app folder
