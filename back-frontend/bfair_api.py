from flask import Flask, request
from flask_sqlalchemy import SQLAlchemy
from flask_migrate import Migrate

app = Flask(__name__)
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

app.config['SQLALCHEMY_DATABASE_URI'] = "postgresql://bfair_user:password@localhost:5432/smart"
db = SQLAlchemy(app)
migrate = Migrate(app, db)


class SampleModel(db.Model):
    __tablename__ = "Sample"

    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    Sample_name = db.Column(db.String())

    def __init__(self, Sample_name):
        self.Sample_name = Sample_name

    def __repr__(self):
        return f"<Sample {self.Sample_name}>"


class MethodologyModel(db.Model):
    __tablename__ = "Methodology"

    id = db.Column(db.Integer, primary_key=True)
    Methodology = db.Column(db.String(120), nullable=False)
    Sample_id = db.Column(db.Integer, db.ForeignKey('Sample.id'),
        nullable=False)