sudo apt update -y && apt upgrade -y
sudo apt-get install python3-pip python3-venv
sudo apt install apache2
sudo apt-get install libapache2-mod-wsgi-py3

mkdir -p /opt/flask
cd /opt/flask/
python3 -m venv venv
source venv/bin/activate
pip install flask
cp start_scripts/app.py /opt/flask/
cp start_scripts/flaskapp.wsgi /opt/flask/
cp start_scripts/flask.conf /etc/apache2/sites-available/flask.conf

sudo a2ensite flask.conf
sudo systemctl restart apache2
