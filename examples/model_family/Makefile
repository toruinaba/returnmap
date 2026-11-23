.PHONY: test check docker-image docker-check docker-run

IMAGE = mymodels-fortran
PACKAGE = mymodels

# Runs as the invoking user, so that what the build writes into the bind mount
# stays writable on the host.
DOCKER_RUN = docker run --rm --user $$(id -u):$$(id -g) \
	-v $$(pwd):/workspace -w /workspace -e HOME=/tmp \
	-e PYTHONPATH=/workspace/src $(IMAGE)

## Tests on this machine.  Without gfortran the UMAT tests are skipped.
test:
	uv run pytest

## The report of every model on this machine.  An unbuilt UMAT reads "not built".
check:
	uv run returnmap check --family $(PACKAGE)

## Build the image (once, and again when pyproject.toml changes).
docker-image:
	docker build -t $(IMAGE) .

## In the container: build every declared UMAT, run the whole test suite, and
## print the report of every model.  An unbuilt UMAT counts as a failure here.
docker-check:
	$(DOCKER_RUN) bash -c "returnmap build --family $(PACKAGE) \
		&& python -m pytest -q -rs \
		&& returnmap check --family $(PACKAGE) --require-built"

## One command in the container:  make docker-run CMD="returnmap check --family mymodels --range"
docker-run:
	$(DOCKER_RUN) bash -c "$(CMD)"
