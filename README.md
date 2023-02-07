# AlphaCharges

**TODO**: Write a project description

## How to run

To run AlphaCharges locally, you will need to have [Python 3](https://www.python.org/downloads/) and [pip](https://pip.pypa.io/en/stable/installing/) installed.

Then, install the project dependencies by running:

```bash
$ sudo python -m venv venv
$ . venv/bin/activate
$ pip install -r requirements.txt
```

**BUG**: Before running the project you will need to copy over this file:

```bash
$ sudo chown -R $USER:$USER /opt
$ mkdir -p /opt/venv/bin
$ cp ./venv/bin/pdb2pqr30 /opt/venv/bin/pdb2pqr
```

Run the project by running the following command inside the virtual environment:

```bash
(venv) $ flask --app app/routes.py run
```

## License
MIT