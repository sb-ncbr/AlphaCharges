from flask import Flask

application = Flask(__name__)

application.jinja_env.trim_blocks = True
application.jinja_env.lstrip_blocks = True

application.config['SECRET_KEY'] = "asdfasdf"
