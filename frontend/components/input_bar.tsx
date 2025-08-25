

"use client";

import type React from "react";
import { useState } from "react";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Textarea } from "@/components/ui/textarea";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Upload, Info, Play, Loader2 } from "lucide-react";
import { ResultsTable, type PredictionResult } from "./ResultsTable";

export default function InputBar() {
  const [activeTab, setActiveTab] = useState("single");
  const [smilesInput, setSmilesInput] = useState("");
  const [multipleSmilesInput, setMultipleSmilesInput] = useState("");
  const [selectedMethod, setSelectedMethod] = useState("");
  const [uploadedFile, setUploadedFile] = useState<File | null>(null);
  const [results, setResults] = useState<PredictionResult[]>([]);
  const [isLoading, setIsLoading] = useState(false);

  const handleSubmit = async () => {
    setIsLoading(true);
    setResults([]);

    const formData = new FormData();
    formData.append("method", selectedMethod);

    let data: string | File | null = null;
    let inputType: string = "";

    if (activeTab === "single") {
      inputType = "smiles";
      data = smilesInput;
    } else {
      if (uploadedFile) {
        inputType = "file";
        data = uploadedFile;
      } else if (multipleSmilesInput.trim()) {
        inputType = "file";
        const smilesArray = multipleSmilesInput.trim().split(/[\s,]+/);
        let csvContent = "compound,smiles\n";
        smilesArray.forEach((s, i) => {
          if (s) csvContent += `Molecule_${i + 1},${s}\n`;
        });
        data = new File([csvContent], "smiles.csv", { type: "text/csv" });
      }
    }

    if (!data) {
      alert("Please provide input.");
      setIsLoading(false);
      return;
    }

    formData.append("inputType", inputType);
    formData.append("data", data);

    try {
      const response = await fetch("/api/predict", {
        method: "POST",
        body: formData,
      });

      const responseData = await response.json();

      if (!response.ok) {
        throw new Error(responseData.error || "An unknown error occurred.");
      }

      setResults(responseData);
    } catch (error) {
      const errorMessage = error instanceof Error ? error.message : "An unexpected network error occurred.";
      setResults([{ smiles: "", error: errorMessage, compound: "Error" }]);
    } finally {
      setIsLoading(false);
    }
  };

  const isSubmitDisabled = !selectedMethod || isLoading || (activeTab === "single" ? !smilesInput : (!uploadedFile && !multipleSmilesInput));

  return (
    <div className="max-w-screen-2xl mx-auto p-6">
      <div className="flex gap-6 h-[calc(100vh-8rem)]">
        {/* Input Section */}
        <Card className="bg-white shadow-sm w-[28%] flex-shrink-0">
          <CardHeader className="pb-4">
            <CardTitle className="flex items-center gap-2 text-lg font-semibold">
              <Upload className="w-5 h-5" />
              Input
            </CardTitle>
          </CardHeader>

          <CardContent className="space-y-6">
            <Tabs value={activeTab} onValueChange={setActiveTab} className="w-full">
              <TabsList className="grid w-full grid-cols-2 bg-gray-100 rounded-lg">
                <TabsTrigger value="single" className="data-[state=active]:bg-white">
                  Single SMILES
                </TabsTrigger>
                <TabsTrigger value="multiple" className="data-[state=active]:bg-white">
                  Multiple SMILES
                </TabsTrigger>
              </TabsList>
              <TabsContent value="single" className="mt-4 space-y-2">
                <Label htmlFor="smiles-input" className="text-sm font-medium">
                  SMILES
                </Label>
                <Textarea
                  id="smiles-input"
                  placeholder="Enter single SMILES string"
                  value={smilesInput}
                  onChange={(e) => setSmilesInput(e.target.value)}
                  className="min-h-[100px] resize-none"
                />
              </TabsContent>
              <TabsContent value="multiple" className="mt-4 space-y-4">
                <div className="border-2 border-dashed border-gray-300 rounded-lg p-6 text-center bg-gray-50 relative">
                  <Upload className="w-10 h-10 text-gray-400 mx-auto mb-3" />
                  <p className="text-sm font-medium text-gray-700">Click to upload or drag CSV file</p>
                  <Input
                    type="file"
                    accept=".csv"
                    onChange={(e) => setUploadedFile(e.target.files?.[0] ?? null)}
                    className="hidden"
                    id="file-upload"
                  />
                  <Label htmlFor="file-upload" className="absolute inset-0 cursor-pointer" />
                </div>
                {uploadedFile && <p className="text-sm text-green-600">File: {uploadedFile.name}</p>}
                <div className="space-y-2">
                  <Label htmlFor="multiple-smiles-input" className="text-sm font-medium">
                    Or enter multiple SMILES (newline separated)
                  </Label>
                  <Textarea
                    id="multiple-smiles-input"
                    placeholder="C1=CC=CC=C1\nCCO"
                    value={multipleSmilesInput}
                    onChange={(e) => setMultipleSmilesInput(e.target.value)}
                    className="min-h-[70px] resize-none"
                  />
                </div>
              </TabsContent>
            </Tabs>

            <div className="space-y-2">
              <div className="flex items-center gap-2">
                <Label htmlFor="method-select" className="text-sm font-medium">
                  Select method
                </Label>
                <Info className="w-4 h-4 text-gray-400" />
              </div>
              <Select value={selectedMethod} onValueChange={setSelectedMethod}>
                <SelectTrigger className="w-full h-10">
                  <SelectValue placeholder="Choose prediction method" />
                </SelectTrigger>
                <SelectContent>
                  <SelectItem value="girolami">Girolami</SelectItem>
                  <SelectItem value="immirzi">Immirzi-Perini</SelectItem>
                </SelectContent>
              </Select>
            </div>

            <Button
              onClick={handleSubmit}
              className="w-full bg-[var(--color-chart-2)] hover:bg-[var(--color-chart-3)] text-white"
              disabled={isSubmitDisabled}
            >
              {isLoading ? <Loader2 className="w-4 h-4 mr-2 animate-spin" /> : <Play className="w-4 h-4 mr-2" />}
              Submit
            </Button>
          </CardContent>
        </Card>

        {/* Results Section */}
        <Card className="bg-white shadow-sm flex-1">
          <CardHeader className="pb-4">
            <CardTitle className="text-lg font-semibold ">
              {selectedMethod && results.length > 0
                ? `${
                    selectedMethod === "immirzi"
                      ? "Immirzi-Perini"
                      : selectedMethod.charAt(0).toUpperCase() + selectedMethod.slice(1)
                  } Method: Results`
                : "Results"}
            </CardTitle>
          </CardHeader>
          <CardContent className="flex items-center justify-center h-[calc(100%-4rem)]">
            <ResultsTable results={results} isLoading={isLoading} />
          </CardContent>
        </Card>
      </div>
    </div>
  );
}