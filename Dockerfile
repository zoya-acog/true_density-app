# Stage 1: Build the Next.js frontend
FROM node:18-slim AS frontend-builder

WORKDIR /app

# Copy frontend package manager files
COPY frontend/package.json frontend/package-lock.json* ./
RUN npm install

# Copy the rest of the frontend code
COPY frontend/ ./

# Build the Next.js application
RUN npm run build

# Stage 2: Setup the Python environment and final application
FROM node:18-slim AS runner

WORKDIR /app

# Install Python and venv
RUN apt-get update && \
    apt-get install -y python3 python3-pip python3-venv && \
    rm -rf /var/lib/apt/lists/*

# Create a virtual environment
RUN python3 -m venv /opt/venv

# Add the venv to the PATH
ENV PATH="/opt/venv/bin:$PATH"

# Copy Python requirements and install dependencies into the venv
COPY True_density/requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# Copy the Python scripts
COPY True_density ./True_density

# Copy the built Next.js application from the builder stage
COPY --from=frontend-builder /app/.next ./.next
COPY --from=frontend-builder /app/public ./public
COPY --from=frontend-builder /app/package.json ./

# Set the command to start the Next.js server
CMD ["npm", "start"]
