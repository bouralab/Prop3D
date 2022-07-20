FROM quay.io/ucsc_cgl/toil:5.5.0

# We make sure to create a work directory; Toil can't hot-deploy a
# script from the root of the filesystem, which is where we start.
RUN mkdir /tmp/work
WORKDIR /tmp/work

# We make a virtual environment to allow workflow dependencies to be
# hot-deployed.
#
# We don't really make use of it in this example, but for workflows
# that depend on PyPI packages we will need this.
#
# We use --system-site-packages so that the Toil installed in the
# appliance image is still available.
ENV VIRTUAL_ENV=/opt/venv
RUN virtualenv --python python3 --system-site-packages $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
RUN python -m pip install awscli botocore==1.23.23
RUN python -m pip install numpy pandas tables scikit-learn Biopython seaborn
RUN python -m pip install tqdm dask dask[dataframe] joblib tornado toolz>=0.7.3
RUN python -m pip install partd>=0.3.8 cloudpickle>=0.2.1 freesasa spython h5py h5pyd
RUN python -m pip install h5py h5pyd
RUN python -m pip install -U urllib3
RUN python -m pip install -U requests


ADD . / /tmp/work/molmimic/
WORKDIR /tmp/work/molmimic/


RUN echo "[easy_install]" > ~/.distutils.cfg
RUN echo "zip_ok = False" >> ~/.distutils.cfg

#Update toil
RUN python3 -m pip uninstall -y toil
WORKDIR /tmp/work/molmimic/toil
RUN python setup.py sdist
RUN pip install /tmp/work/molmimic/toil/dist/toil-5.6.0a1.tar.gz[all]

WORKDIR /tmp/work/molmimic/

#Now we will install molmimic into the virtualenv
RUN rm requirements.txt
RUN touch requirements.txt
RUN python -m pip install .
#RUN pip install /tmp/work/molmimic/dist/molmimic-0.0.1.tar.gz
WORKDIR /tmp/work
#RUN rm -rf molmimic

RUN mkdir /tmp/scratch

ENTRYPOINT ["python", "-m", "molmimic.generate_data.main"]
