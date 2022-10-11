#  Source: ubuntu-jammy-x86_64
#  Flavor: standard.large
#  Networks: 147-251-115-pers-proj-net

# https://techtutorguro.com/how-to-install-flask-on-ubuntu-22-04-with-apache/



#sudo apt update
#sudo apt install git
#sudo bash start.sh
# git clone https://github.com/dargen3/AlphaCharges in /opt folder



sudo apt update -y && apt upgrade -y
sudo apt-get install python3-pip python3-venv
sudo apt install apache2
sudo apt-get install libapache2-mod-wsgi-py3

python3 -m venv venv
source venv/bin/activate
pip install flask requests scipy rdkit numba numpy pdb2pqr
sudo rm -f /etc/apache2/sites-available/*
cp /home/ubuntu/AlphaCharges/start_scripts/alpha_charges.conf /etc/apache2/sites-available/
sudo chown -R www-data:www-data /opt
sudo chmod o+rx app/AlphaCharges.wsgi
sudo chmod o+rx app/routes.py



sudo a2ensite AlphaCharges.conf
sudo systemctl restart apache2

# cat /var/log/apache2/error.log0
