<div align="center">
  <br>
  <br>
  <a href="https://github.com/sb-ncbr/AlphaCharges"><img src="https://github.com/sb-ncbr/AlphaCharges/blob/50265b26f8748e4afa3b9d4619e8f04e83640b13/app/static/assets/logo.png" alt="AlphaCharges" width="220"></a>
  <br>
  <br>
</div>

[αCharges](https://alphacharges.ncbr.muni.cz/) is a web application for the calculation of partial atomic charges on protein structures predicted by the [AlphaFold2](https://www.nature.com/articles/s41586-021-03819-2) algorithm and deposited in the [AlphaFoldDB](https://academic.oup.com/nar/article/50/D1/D439/6430488) database. The charges are computed by the [SQE+qp](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00528-w) empirical method, which quality is comparable to quantum mechanical charge calculation approaches. Before computation of the charges, αCharges protonates the input protein structures by [PROPKA3](https://pubs.acs.org/doi/full/10.1021/ct100578z). The details about the methodology and usage are described in the [manual](https://github.com/sb-ncbr/AlphaCharges/wiki). This website is free and open to all users and there is no login requirement.

## How to run

To run AlphaCharges locally, you will need to have [Python 3.9](https://www.python.org/downloads/) and [pip](https://pip.pypa.io/en/stable/installing/) installed.

Then, clone project and install the project dependencies by running:

```bash
$ cd /opt
$ git clone https://github.com/sb-ncbr/AlphaCharges
$ sudo python3.9 -m venv venv
$ . venv/bin/activate
$ pip install -r requirements.txt
```
Run the project by running the following command inside the virtual environment:

```bash
(venv) $ cd /opt/Alphacharges/app
(venv) $ export FLASK_APP=routes.py
(venv) $ flask run
```
and point your browser to localhost:5000/.

## How to cite

If you found AlphaCharges helpful, please cite: [Schindler, O., Berka, K., Cantara, A., Křenek, A., Tichý, D., Raček, T., & Svobodová, R. (2023). αCharges: Partial atomic charges for AlphaFold structures in high quality. Nucleic Acids Research.](https://doi.org/10.1093/nar/gkad349)

## License
MIT
