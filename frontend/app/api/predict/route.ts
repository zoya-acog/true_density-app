import { NextRequest, NextResponse } from 'next/server';
import { spawn } from 'child_process';
import path from 'path';
import fs from 'fs/promises';
import os from 'os';

// Define the structure of the prediction results
interface PredictionResult {
  compound?: string; // Optional for single SMILES
  smiles: string;
  density?: number;
  error_density?: number;
  volume?: number;
  image?: string; // base64 encoded image
  error?: string;
}

// Path to the Python executable, 'python3' should be in the system's PATH in the Docker container
const PYTHON_EXECUTABLE = 'python3';

// Path to the CLI script, resolved from the current working directory
const CLI_SCRIPT_PATH = path.resolve(process.cwd(), 'True_density/cli.py');

async function runCliProcess(args: string[]): Promise<PredictionResult[]> {
  return new Promise((resolve, reject) => {
    const commandArgs = [CLI_SCRIPT_PATH, ...args, '--json'];
    const pythonProcess = spawn(PYTHON_EXECUTABLE, commandArgs);

    let stdout = '';
    let stderr = '';

    pythonProcess.stdout.on('data', (data) => {
      stdout += data.toString();
    });

    pythonProcess.stderr.on('data', (data) => {
      stderr += data.toString();
    });

    pythonProcess.on('close', (code) => {
      if (code !== 0) {
        console.error(`Python script exited with code ${code}`);
        console.error(`Stderr: ${stderr}`);
        // Try to parse stderr for a JSON error message from the script
        try {
            const errorJson = JSON.parse(stderr);
            return reject(new Error(errorJson.error || 'An unknown error occurred in the Python script.'));
        } catch (e) {
            return reject(new Error(stderr || 'An unknown error occurred in the Python script.'));
        }
      }

      try {
        const results: PredictionResult[] = JSON.parse(stdout);
        resolve(results);
      } catch (error) {
        console.error('Failed to parse JSON from python script:', stdout);
        reject(new Error('Failed to parse response from the prediction script.'));
      }
    });

    pythonProcess.on('error', (error) => {
        console.error('Failed to start subprocess.', error);
        reject(new Error('Failed to start the prediction process.'));
    });
  });
}

export async function POST(req: NextRequest) {
  try {
    const formData = await req.formData();
    const method = formData.get('method') as string; // e.g., 'girolami' or 'immirzi'
    const inputType = formData.get('inputType') as string; // e.g., 'smiles', 'name', or 'file'
    const data = formData.get('data');

    if (!method || !inputType || !data) {
      return NextResponse.json({ error: 'Missing required fields: method, inputType, data' }, { status: 400 });
    }

    let cliArgs: string[];

    if (inputType === 'file') {
      const file = data as File;
      const tempFilePath = path.join(os.tmpdir(), file.name);
      const fileBuffer = Buffer.from(await file.arrayBuffer());
      await fs.writeFile(tempFilePath, fileBuffer);
      
      cliArgs = [`${method}-file`, tempFilePath];
    } else {
      const textData = data as string;
      cliArgs = [`${method}-${inputType}`, textData];
    }

    const results = await runCliProcess(cliArgs);

    // Clean up temp file if one was created
    if (inputType === 'file') {
        const file = data as File;
        const tempFilePath = path.join(os.tmpdir(), file.name);
        await fs.unlink(tempFilePath);
    }

    // Ensure all fields are present, excluding compound for single SMILES
    const processedResults = results.map((item: any, index) => {
      const result: PredictionResult = {
        smiles: item.smiles || 'N/A',
        image: item.image,
        density: item.density,
        error_density: item.error_density,
        volume: item.volume,
        error: item.error,
      };
      if (inputType !== 'smiles') {
        result.compound = item.compound || `Molecule_${index + 1}`;
      }
      return result;
    });

    return NextResponse.json(processedResults);

  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'An unexpected error occurred.';
    return NextResponse.json({ error: errorMessage }, { status: 500 });
  }
}