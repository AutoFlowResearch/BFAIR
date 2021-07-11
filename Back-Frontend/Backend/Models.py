from flask_sqlalchemy import SQLAlchemy
from flask_migrate import Migrate
from flask import Flask, jsonify, request
from flask_marshmallow import Marshmallow, Schema
from flask_cors import CORS

"""Database models."""
# from . import db
from flask_login import UserMixin, LoginManager
from werkzeug.security import generate_password_hash, check_password_hash

app = Flask(__name__)
ma = Marshmallow(app)
CORS(app)
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

app.config['SQLALCHEMY_DATABASE_URI'] = "postgresql://bfair_user:password@localhost:5432/smart"
login = LoginManager()
db = SQLAlchemy(app)
migrate = Migrate(app, db)
db.init_app(app)
login.init_app(app)
login.login_view = 'login'


class UserModel(UserMixin, db.Model):
    """User account model."""

    __tablename__ = 'users'
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(100), nullable=False, unique=False)
    email = db.Column(db.String(40), unique=True, nullable=False)
    password = db.Column(db.String(200), primary_key=False, unique=False, nullable=False)

    # created_on = db.Column(db.DateTime, index=False, unique=False, nullable=True)
    # last_login = db.Column(db.DateTime, index=False, unique=False, nullable=True)
    def __repr__(self):
        return f"<Sample {self.name}>"

    def set_password(self, password):
        """Create hashed password."""
        self.password = generate_password_hash(
            password,
            method='sha256'
        )

    def check_password(self, password):
        """Check hashed password."""
        return check_password_hash(self.password, password)

    def __repr__(self):
        return '<User {}>'.format(self.username)

# @login.user_loader
# def load_user(id):
#     return UserModel.query.get(int(id))
