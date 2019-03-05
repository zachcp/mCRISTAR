"""`main` is the top level module for your Flask application."""

# [START gae_python37_app]

# Import the Flask Framework
from io import StringIO
from io import BytesIO
import json

from Bio import SeqIO
from flask import Flask, render_template, request

from mCRISTAR.data import selective_promoters, nonselective_promoters
from mCRISTAR.mCRISTAR import GBKProcessor, GapProcessor

app = Flask(__name__)

@app.route('/')
def hello():
    return render_template('home.html')

@app.route('/data')
def data():
    return render_template('data.html',
                           selective_promoters=selective_promoters,
                           nonselective_promoters=nonselective_promoters)

@app.route("/refactor", methods=["GET"])
def refactor():
    return render_template('refactor.html')

@app.route("/upload", methods=["GET","POST"])
def upload():
    "upload and process the GBK file"
    upload_file  = request.files["file"]
    genbank_string  = upload_file.read()
    gb_stringhandle = StringIO(str(genbank_string, 'utf-8'))
    gbk = SeqIO.read(gb_stringhandle,"gb")
    gbkProcessor = GBKProcessor(gbk=gbk)
    return json.dumps(gbkProcessor.export())

@app.route("/makecassettes", methods=["GET","POST"])
def makecassettes():
    "find crisprsites and create primers"
    gapdata = request.json['gaps']
    cassetes = GapProcessor(gapdata)
    return json.dumps(cassetes.export())

@app.errorhandler(404)
def page_not_found(e):
    """Return a custom 404 error."""
    return 'Sorry, Nothing at this URL.', 404

@app.errorhandler(500)
def application_error(e):
    """Return a custom 500 error."""
    return 'Sorry, unexpected error: {}'.format(e), 500

if __name__ == "__main__":
    app.debug=True
    app.run()

# [END gae_python37_app]
