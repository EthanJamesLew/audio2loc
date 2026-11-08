# audio2loc: Doppler/ADS-B geolocation CLI + MCP server.
#
#   docker build -t audio2loc .                      # slim image (XML .qgs projects)
#   docker build -t audio2loc:qgis --target qgis .   # + native QGIS API (.qgz export)
#   docker build -t audio2loc:cuda --target cuda .   # + CUDA kernels (run with --device nvidia.com/gpu=all)
#
#   docker run --rm -i -v audio2loc-db:/data -v "$PWD":/work audio2loc            # MCP server on stdio
#   docker run --rm -v audio2loc-db:/data audio2loc update --latest --bbox -125,24,-66,50
#   docker run --rm -v audio2loc-db:/data -v "$PWD":/work audio2loc locate --audio /work/clip.mp4 ...

ARG UBUNTU=24.04
ARG CUDA_IMAGE=12.6.3

# ---------------------------------------------------------------- build ---
FROM ubuntu:${UBUNTU} AS build
RUN apt-get update && apt-get install -y --no-install-recommends g++ cmake make zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /src
ARG AUDIO2LOC_VERSION=
COPY CMakeLists.txt .
COPY cmake cmake
COPY src src
# portable binary: no -march=native inside images
RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DAUDIO2LOC_NATIVE=OFF -DAUDIO2LOC_BUILD_TESTS=OFF -DAUDIO2LOC_VERSION="$AUDIO2LOC_VERSION" && cmake --build build -j"$(nproc)"

# -------------------------------------------------------------- runtime ---
FROM ubuntu:${UBUNTU} AS runtime
LABEL org.opencontainers.image.title="audio2loc" \
      org.opencontainers.image.description="Geolocate a recording from aircraft Doppler curves and ADS-B history"
ARG YTDLP_VERSION=latest
RUN apt-get update && apt-get install -y --no-install-recommends \
        ffmpeg curl ca-certificates tar python3 zlib1g \
    && rm -rf /var/lib/apt/lists/* \
    && if [ "$YTDLP_VERSION" = latest ]; then u=https://github.com/yt-dlp/yt-dlp/releases/latest/download/yt-dlp; \
       else u=https://github.com/yt-dlp/yt-dlp/releases/download/$YTDLP_VERSION/yt-dlp; fi \
    && curl -fsSL "$u" -o /usr/local/bin/yt-dlp && chmod +x /usr/local/bin/yt-dlp
COPY --from=build /src/build/audio2loc /usr/local/bin/audio2loc
# ubuntu:24.04 ships user "ubuntu" (uid 1000), matching typical host bind-mount ownership
RUN mkdir -p /data /work && chown ubuntu:ubuntu /data /work
USER ubuntu
ENV AUDIO2LOC_DB=/data AUDIO2LOC_WORK=/work
VOLUME ["/data", "/work"]
WORKDIR /work
ENTRYPOINT ["audio2loc"]
CMD ["mcp"]

# ------------------------------------------------- QGIS-enabled variant ---
FROM ubuntu:${UBUNTU} AS build-qgis
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        g++ cmake make zlib1g-dev libqgis-dev qtbase5-dev libqt5xml5 \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /src
ARG AUDIO2LOC_VERSION=
COPY CMakeLists.txt .
COPY cmake cmake
COPY src src
RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DAUDIO2LOC_NATIVE=OFF -DAUDIO2LOC_BUILD_TESTS=OFF -DAUDIO2LOC_VERSION="$AUDIO2LOC_VERSION" -DAUDIO2LOC_WITH_QGIS=ON \
    && cmake --build build -j"$(nproc)"

FROM runtime AS qgis
USER root
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends libqgis-core3.34.4 \
    || (apt-get install -y --no-install-recommends $(apt-cache search --names-only '^libqgis-core3' | cut -d' ' -f1 | head -1)) \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build-qgis /src/build/audio2loc /usr/local/bin/audio2loc
ENV QT_QPA_PLATFORM=offscreen
USER ubuntu

# ------------------------------------------------- CUDA-enabled variant ---
#   docker build -t audio2loc:cuda --target cuda .
#   docker run --rm -i --device nvidia.com/gpu=all -v audio2loc-db:/data -v "$PWD":/work audio2loc:cuda
#   (or --gpus all on hosts using the legacy nvidia runtime). Without a GPU it falls back to CPU.
FROM nvidia/cuda:${CUDA_IMAGE}-devel-ubuntu24.04 AS build-cuda
RUN apt-get update && apt-get install -y --no-install-recommends g++ cmake make zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /src
ARG AUDIO2LOC_VERSION=
COPY CMakeLists.txt .
COPY cmake cmake
COPY src src
# Turing .. Hopper SASS + PTX for newer GPUs
RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DAUDIO2LOC_NATIVE=OFF -DAUDIO2LOC_BUILD_TESTS=OFF -DAUDIO2LOC_VERSION="$AUDIO2LOC_VERSION" -DAUDIO2LOC_CUDA=ON \
        -DCMAKE_CUDA_ARCHITECTURES="75-real;80-real;86-real;89-real;90" \
    && cmake --build build -j"$(nproc)"

FROM nvidia/cuda:${CUDA_IMAGE}-runtime-ubuntu24.04 AS cuda
ARG YTDLP_VERSION=latest
RUN apt-get update && apt-get install -y --no-install-recommends \
        ffmpeg curl ca-certificates tar python3 zlib1g \
    && rm -rf /var/lib/apt/lists/* \
    && if [ "$YTDLP_VERSION" = latest ]; then u=https://github.com/yt-dlp/yt-dlp/releases/latest/download/yt-dlp; \
       else u=https://github.com/yt-dlp/yt-dlp/releases/download/$YTDLP_VERSION/yt-dlp; fi \
    && curl -fsSL "$u" -o /usr/local/bin/yt-dlp && chmod +x /usr/local/bin/yt-dlp
COPY --from=build-cuda /src/build/audio2loc /usr/local/bin/audio2loc
RUN mkdir -p /data /work && chown ubuntu:ubuntu /data /work
USER ubuntu
ENV AUDIO2LOC_DB=/data AUDIO2LOC_WORK=/work NVIDIA_DRIVER_CAPABILITIES=compute,utility
VOLUME ["/data", "/work"]
WORKDIR /work
ENTRYPOINT ["audio2loc"]
CMD ["mcp"]

# default target (last stage): the slim runtime
FROM runtime AS slim
