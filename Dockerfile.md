# Dockerfile for Glimepiride Model

## Build image
To build the latest image use:
```bash
docker build -f Dockerfile -t matthiaskoenig/dapagliflozin:0.9.8 -t matthiaskoenig/dapagliflozin:latest .
```

## Push images
The image is pushed to dockerhub: [Docker Hub – Dapagliflozin](https://hub.docker.com/repository/docker/matthiaskoenig/dapagliflozin/general)

```
docker login
docker push --all-tags matthiaskoenig/dapagliflozin
```

## Run container
To use the latest container version interactively use:

```bash
docker run -v "${PWD}/results:/results" -it matthiaskoenig/dapagliflozin:latest /bin/bash
```

To use a specific container version provide the version tag:
```bash
docker run -v "${PWD}/results:/results" -it matthiaskoenig/dapagliflozin:0.9.8 /bin/bash
```

## Run simulations
Run the complete analysis:
```bash
uv run run_dapagliflozin -a all -r /results
```
The results are written into the mounted `/results` folder on the host.

In case of permission issues with the mounted folder, adjust ownership and access rights with:
```bash
sudo chown $(id -u):$(id -g) -R "${PWD}/results"
sudo chmod 775 "${PWD}/results"
```
