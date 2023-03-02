FROM globalseismology/avni-buildenv-jammy:latest

EXPOSE 8888

RUN useradd \
  --create-home \
  avni_user

USER avni_user

WORKDIR /home/avni_user

RUN git clone 'https://github.com/globalseismology/avni.git'
ENV AVNI_DIR /home/avni_user/avni
WORKDIR /home/avni_user/avni

RUN cd $AVNI_DIR && pip3 install .
