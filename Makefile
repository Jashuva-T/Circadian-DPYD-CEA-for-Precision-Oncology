.PHONY: install test run clean
install:
	pip install -r requirements.txt && pip install -e .
test:
	pytest -q
run:
	python -m circadian_dpyd.pipeline --config config/params.yaml --out results
clean:
	rm -rf results/figures/* results/tables/* results/summary.json
