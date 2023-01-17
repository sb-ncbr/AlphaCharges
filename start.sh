# sudo apt update -y && sudo apt upgrade -y
# sudo apt install git
# sudo git clone https://github.com/dargen3/AlphaCharges in /opt folder




sudo apt-get install python3-pip python3-venv apache2 libapache2-mod-wsgi-py3 gemmi

sudo python3 -m venv venv
source venv/bin/activate
sudo chown -R ubuntu:ubuntu /opt
pip install flask requests scipy numba numpy pdb2pqr scikit-learn rdkit==2022.03.5
sudo rm -f /etc/apache2/sites-available/*
sudo cp AlphaCharges.conf /etc/apache2/sites-available/
sudo chown -R www-data:www-data /opt
sudo chmod o+rx app/AlphaCharges.wsgi
sudo chmod o+rx app/routes.py


sudo a2ensite AlphaCharges.conf

sudo a2enmod ssl
sudo a2enmod brotli
sudo a2enmod http2

sudo systemctl restart apache2

# after update
# sudo chown -R www-data:www-data /opt
# cat /var/log/apache2/error.log

