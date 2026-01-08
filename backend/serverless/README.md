# RunPod Serverless Deployment Guide

Deploy the Foundry Protein Design backend as a serverless endpoint for **90%+ cost savings** compared to always-on GPU pods.

## Architecture

```
Frontend (Vercel) ─── Edge Function ───> RunPod Serverless ───> Network Volume
                      /api/runpod/*       (Flex Workers)        (Checkpoints)
                           │
                           └───> Supabase (Job History)
```

## Cost Comparison

| Setup | Monthly Cost |
|-------|-------------|
| Always-On Pod (A40) | ~$570/month |
| **Serverless (10 hrs/week)** | **~$27/month** |

## Prerequisites

1. **RunPod Account** with API access
2. **Docker Hub** account (for pushing images)
3. **Supabase** account (free tier for job persistence)
4. Model checkpoints (~20GB)

## Step 1: Build and Push Docker Image

```bash
cd backend/serverless

# Build the image
docker build -t your-username/foundry-serverless:latest .

# Push to Docker Hub
docker push your-username/foundry-serverless:latest
```

> **Note:** The first build takes ~30 minutes due to Foundry installation.

## Step 2: Create Network Volume

1. Go to [RunPod Console → Storage](https://www.runpod.io/console/storage)
2. Click **New Network Volume**
3. Configure:
   - **Name:** `foundry-checkpoints`
   - **Size:** 100GB
   - **Region:** Choose closest to your users (e.g., `US-OR-1`)
4. Click **Create**

### Upload Checkpoints

Use RunPod's S3-compatible API:

```bash
# Get credentials from RunPod Console → Settings → API Keys → S3 Credentials
export AWS_ACCESS_KEY_ID=<your_access_key>
export AWS_SECRET_ACCESS_KEY=<your_secret_key>
export S3_ENDPOINT=https://us-or-1.storage.runpod.io  # Match your region

# Upload checkpoints (from a machine with the models)
aws s3 sync /workspace/checkpoints s3://foundry-checkpoints/checkpoints/ \
    --endpoint-url $S3_ENDPOINT

# Verify upload
aws s3 ls s3://foundry-checkpoints/checkpoints/ --endpoint-url $S3_ENDPOINT
```

**Expected structure:**
```
/runpod-volume/checkpoints/
├── rfd3/
├── rf3/
└── mpnn/
```

## Step 3: Create Serverless Endpoint

1. Go to [RunPod Console → Serverless](https://www.runpod.io/console/serverless)
2. Click **New Endpoint**
3. Configure:

| Setting | Value |
|---------|-------|
| **Name** | `foundry-protein-design` |
| **Template Source** | Docker Image |
| **Docker Image** | `your-username/foundry-serverless:latest` |
| **GPU** | A40 (48GB) or L40 |
| **Min Workers** | 0 (scale to zero) |
| **Max Workers** | 2 |
| **Idle Timeout** | 30 seconds |
| **Execution Timeout** | 3600 seconds |
| **Network Volume** | `foundry-checkpoints` |

4. Click **Create Endpoint**
5. Copy your **Endpoint ID** (e.g., `abc123xyz`)

## Step 4: Get API Key

1. Go to [RunPod Console → Settings → API Keys](https://www.runpod.io/console/settings/api-keys)
2. Create a new key with **Serverless** permissions
3. Copy the API key (e.g., `rp_XXXX...`)

## Step 5: Test the Endpoint

```bash
# Test health check
curl -X POST https://api.runpod.ai/v2/{ENDPOINT_ID}/run \
  -H "Authorization: Bearer {API_KEY}" \
  -H "Content-Type: application/json" \
  -d '{"input": {"task": "health"}}'

# Get result (replace JOB_ID from response)
curl https://api.runpod.ai/v2/{ENDPOINT_ID}/status/{JOB_ID} \
  -H "Authorization: Bearer {API_KEY}"
```

**Expected first request:** ~60-90 seconds (cold start)
**Subsequent requests:** ~1-5 seconds

## Step 6: Configure Frontend

Add environment variables to your Vercel project:

```env
RUNPOD_API_KEY=rp_XXXX...
RUNPOD_ENDPOINT_ID=abc123xyz
```

## API Reference

### Request Format

```json
{
  "input": {
    "task": "rfd3" | "rf3" | "mpnn" | "rmsd" | "health",
    // Task-specific parameters below
  }
}
```

### RFD3 Design

```json
{
  "input": {
    "task": "rfd3",
    "contig": "100",
    "num_designs": 1,
    "seed": 42
  }
}
```

### RF3 Prediction

```json
{
  "input": {
    "task": "rf3",
    "sequence": "MSKGEELFT...",
    "name": "my_protein"
  }
}
```

### MPNN Sequence Design

```json
{
  "input": {
    "task": "mpnn",
    "pdb_content": "ATOM...",
    "num_sequences": 8,
    "temperature": 0.1,
    "model_type": "ligand_mpnn"
  }
}
```

### RMSD Validation

```json
{
  "input": {
    "task": "rmsd",
    "pdb_content_1": "ATOM...",
    "pdb_content_2": "ATOM...",
    "backbone_only": true
  }
}
```

### Response Format

```json
{
  "id": "job-123",
  "status": "COMPLETED",
  "output": {
    "status": "completed",
    "result": { ... }
  }
}
```

**Status values:** `IN_QUEUE`, `IN_PROGRESS`, `COMPLETED`, `FAILED`

## Troubleshooting

### Cold Start Too Long (>2 min)

- Check if network volume is properly mounted
- Verify checkpoints are in correct location
- Consider keeping 1 active worker ($150/month)

### "Module not found" Errors

- Rebuild Docker image with latest Foundry
- Check FOUNDRY_CHECKPOINT_DIRS environment variable

### GPU Out of Memory

- Switch to A100 (80GB) for large proteins
- Reduce batch size / num_designs

### Job Stuck in IN_QUEUE

- Check worker logs in RunPod console
- Verify GPU availability in your region

## Local Testing

```bash
cd backend/serverless

# Install dependencies
pip install runpod

# Run local server (mimics RunPod)
python handler.py --rp_serve_api

# Test locally
curl -X POST http://localhost:8000/run \
  -H "Content-Type: application/json" \
  -d '{"input": {"task": "health"}}'
```

## Monitoring

- **RunPod Dashboard:** Real-time worker status, logs, metrics
- **Supabase Dashboard:** Job history and persistence
- **Vercel Dashboard:** Edge function invocations

## Cost Optimization Tips

1. **Scale to zero:** Set min workers = 0
2. **Short idle timeout:** 30 seconds is optimal
3. **Use A40 over A100:** Sufficient for most proteins, 50% cheaper
4. **Batch requests:** Submit multiple designs in one call
