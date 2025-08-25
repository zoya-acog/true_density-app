# Stage 1: Build the frontend
FROM node:18-alpine AS builder
WORKDIR /app/frontend
COPY frontend/package*.json ./
RUN npm install
COPY frontend/ .
RUN npm run build

# Stage 2: Set up the Python backend
FROM continuumio/miniconda3 AS python-backend
WORKDIR /app
COPY True_density/ ./True_density/
RUN conda install -c conda-forge --yes rdkit &&     pip install -r True_density/requirements.txt

# Stage 3: Final image
FROM node:18-alpine
WORKDIR /app
COPY --from=builder /app/frontend/.next ./frontend/.next
COPY --from=builder /app/frontend/public ./frontend/public
COPY --from=builder /app/frontend/package.json ./frontend/
COPY --from=python-backend /app/True_density/ ./True_density/
COPY --from=python-backend /opt/conda/ /opt/conda/

ENV PATH /opt/conda/bin:$PATH

EXPOSE 3000

CMD ["npm", "start", "--", "--hostname", "0.0.0.0"]

