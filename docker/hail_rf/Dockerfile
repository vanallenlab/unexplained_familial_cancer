# The Germline Genomics of Cancer (G2C)
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Simple Dockerfile for vanallenlab/hail_rf


FROM hailgenetics/hail:0.2.127.post1-py3.11
MAINTAINER "Noah Fields <Noah_Fields@dfci.harvard.edu>"


COPY apply_rf_model.py /opt/

# Launch bash at runtime
CMD ["/bin/bash"]
