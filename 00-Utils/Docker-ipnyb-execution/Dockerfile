# Use an official Python runtime as a parent image
FROM python:3.11

# Set the working directory in the container
WORKDIR /app

# Copy requirements file and script into the container
COPY /00-Utils/Docker-ipnyb-execution/requirements.txt /app/requirements.txt

# Install dependencies
RUN pip install -r requirements.txt

