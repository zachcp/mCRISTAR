"""`main` is the top level module for your Flask application."""

# Import the Flask Framework
import StringIO
import json

from Bio import SeqIO
from flask import Flask, render_template, request

from mCRISTAR.data import selective_promoters, nonselective_promoters
from mCRISTAR.mCRISTAR import mCRISTAR, CassetteCheck, find_crispr_site_JSON, create_promoter_primers_from_JSON, make_crispr_cassette, get_gaps, processGBK, processGaps, processCrisprCassettes

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
    upload_file = request.files["file"]
    genbank_string = upload_file.read()
    gb_stringhandle = StringIO.StringIO(genbank_string)
    gbk = SeqIO.read(gb_stringhandle,"gb")
    gaps = get_gaps(gbk=gbk, mindist=50)
    gaps = processGaps(gaps=gaps, gbk=gbk)
    genes = processGBK(gbk)
    return json.dumps({"genes":genes,
                       "gaps":gaps})

@app.route("/makecassettes", methods=["GET","POST"])
def makecassettes():
    gapdata = request.json['gaps']

    #find crisprsites and create primers
    cassetes = CassetteCheck(gapdata)
    primers = create_promoter_primers_from_JSON(gapdata, selective_promoters, nonselective_promoters, overlaplegth=40)

    return json.dumps({"casscass_dataass_data,
                       "primers":primers})


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