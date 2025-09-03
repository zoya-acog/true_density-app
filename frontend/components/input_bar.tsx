// "use client";

// import type React from "react";
// import { useState, useEffect } from "react";
// import { Button } from "@/components/ui/button";
// import { Input } from "@/components/ui/input";
// import { Label } from "@/components/ui/label";
// import { Textarea } from "@/components/ui/textarea";
// import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
// import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
// import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
// import { Upload, Info, Play, Loader2 } from "lucide-react";
// import { ResultsTable, type PredictionResult } from "./ResultsTable";

// export default function InputBar() {
//   const [activeTab, setActiveTab] = useState("single");
//   const [smilesInput, setSmilesInput] = useState("");
//   const [multipleSmilesInput, setMultipleSmilesInput] = useState("");
//   const [selectedMethod, setSelectedMethod] = useState("");
//   const [uploadedFile, setUploadedFile] = useState<File | null>(null);
//   const [results, setResults] = useState<PredictionResult[]>([]);
//   const [isLoading, setIsLoading] = useState(false);

//   // Effect to clear results when switching tabs
//   useEffect(() => {
//     setResults([]); // Clear results whenever activeTab changes
//   }, [activeTab]);

//   const handleSubmit = async () => {
//     setIsLoading(true);
//     setResults([]);

//     const formData = new FormData();
//     formData.append("method", selectedMethod);

//     let data: string | File | null = null;
//     let inputType: string = "";

//     if (activeTab === "single") {
//       inputType = "smiles";
//       data = smilesInput;
//     } else {
//       if (uploadedFile) {
//         inputType = "file";
//         data = uploadedFile;
//       } else if (multipleSmilesInput.trim()) {
//         inputType = "file";
//         const smilesArray = multipleSmilesInput.trim().split('\n');
//         let csvContent = "compound,smiles\n";
//         smilesArray.forEach(line => {
//           const [compound, smiles] = line.split(',').map(s => s.trim());
//           if (smiles) {
//             csvContent += `${compound || `Molecule_${smilesArray.indexOf(line) + 1}`},${smiles}\n`;
//           }
//         });
//         data = new File([csvContent], "smiles.csv", { type: "text/csv" });
//       }
//     }

//     if (!data) {
//       alert("Please provide input.");
//       setIsLoading(false);
//       return;
//     }

//     formData.append("inputType", inputType);
//     formData.append("data", data);

//     try {
//       const response = await fetch("/api/predict", {
//         method: "POST",
//         body: formData,
//       });

//       const responseData = await response.json();

//       if (!response.ok) {
//         throw new Error(responseData.error || "An unknown error occurred.");
//       }

//       setResults(responseData);
//     } catch (error) {
//       const errorMessage = error instanceof Error ? error.message : "An unexpected network error occurred.";
//       setResults([{ smiles: "", error: errorMessage }]); // Exclude compound
//     } finally {
//       setIsLoading(false);
//     }
//   };

//   const isSubmitDisabled = !selectedMethod || isLoading || (activeTab === "single" ? !smilesInput : (!uploadedFile && !multipleSmilesInput));

//   const handleTextareaFocus = () => {
//     setUploadedFile(null); // Clear uploaded file when focusing on textarea
//   };

//   const handleFileInputChange = (e: React.ChangeEvent<HTMLInputElement>) => {
//     setMultipleSmilesInput(""); // Clear textarea when uploading a file
//     setUploadedFile(e.target.files?.[0] ?? null);
//   };

//   return (
//     <div className="max-w-screen-2xl mx-auto p-6">
//       <div className="flex gap-6 h-[calc(100vh-8rem)]">
//         {/* Input Section */}
//         <Card className="bg-white shadow-sm w-[28%] flex-shrink-0">
//           <CardHeader className="pb-4">
//             <CardTitle className="flex items-center gap-2 text-lg font-semibold">
//               <Upload className="w-5 h-5" />
//               Input
//             </CardTitle>
//           </CardHeader>

//           <CardContent className="space-y-6">
//             <Tabs value={activeTab} onValueChange={setActiveTab} className="w-full">
//               <TabsList className="grid w-full grid-cols-2 bg-gray-100 rounded-lg">
//                 <TabsTrigger value="single" className="data-[state=active]:bg-white">
//                   Single SMILES
//                 </TabsTrigger>
//                 <TabsTrigger value="multiple" className="data-[state=active]:bg-white">
//                   Multiple SMILES
//                 </TabsTrigger>
//               </TabsList>
//               <TabsContent value="single" className="mt-4 space-y-2">
//                 <Label htmlFor="smiles-input" className="text-sm font-medium">
//                   SMILES
//                 </Label>
//                 <Textarea
//                   id="smiles-input"
//                   placeholder="Enter single SMILES string"
//                   value={smilesInput}
//                   onChange={(e) => setSmilesInput(e.target.value)}
//                   className="min-h-[100px] resize-none"
//                 />
//               </TabsContent>
//               <TabsContent value="multiple" className="mt-4 space-y-4">
//                 <div className="space-y-2">
//                   <div className="flex items-center gap-2">
//                     <Label className="text-sm font-medium">Upload csv</Label>
//                     <div className="group relative">
//                       <Info className="w-4 h-4 text-gray-400 cursor-help" />
//                       <span className="absolute hidden group-hover:block bg-gray-500 text-white text-xs rounded p-2 w-64 -top-2 left-6 z-10">
//                         The compound IDs could be their generic names or any other custom identifiers given by the user.
//                       </span>
//                     </div>
//                   </div>
//                   <div className="border-2 border-dashed border-gray-300 rounded-lg p-6 text-center bg-gray-50 relative">
//                     <Upload className="w-10 h-10 text-gray-400 mx-auto mb-3 " />
//                     <p className="text-sm font-medium text-gray-700">Click to upload or drag CSV file with multiple SMILES along with their compound IDs (new line separated).

// </p>
//                     <Input
//                       type="file"
//                       accept=".csv"
//                       onChange={handleFileInputChange}
//                       className="hidden"
//                       id="file-upload"
//                     />
//                     <Label htmlFor="file-upload" className="absolute inset-0 cursor-pointer" />
//                   </div>
//                   {uploadedFile && <p className="text-sm text-green-600">File: {uploadedFile.name}</p>}
//                 </div>
//                 <div className="space-y-2">
//                   <div className="flex items-center gap-2">
//                     <Label htmlFor="multiple-smiles-input" className="text-sm font-medium">
//                       Or enter multiple SMILES (newline separated)
//                     </Label>
//                     <div className="group relative">
//                       <Info className="w-4 h-4 text-gray-400 cursor-help" />
//                       <span className="absolute hidden group-hover:block bg-gray-500 text-white text-xs rounded p-2 w-64 -top-2 left-6 z-10">
//                         The compound IDs could be their generic names or any other custom identifiers given by the user.
//                       </span>
//                     </div>
//                   </div>
//                   <Textarea
//                     id="multiple-smiles-input"
//                     placeholder={`Phenacetin, CCOc1ccc(NC(C)=O)cc1
// AGAN_999, Cc1cc(ccc1OCCCCCCN(S(=O)(=O)c1ccc(cc1)C(C)(C)C)C)S(=O)(=O)N(C)C
// 1, COc1ccc2c(c1)c(c(c(=O)n2C)CC(=O)O)Cl
// A, CC(C(=O)O)c1ccc(cc1)C(=O)C`}
//                     value={multipleSmilesInput}
//                     onChange={(e) => setMultipleSmilesInput(e.target.value)}
//                     onFocus={handleTextareaFocus}
//                     className="min-h-[70px] resize-none"
//                   />
//                 </div>
//               </TabsContent>
//             </Tabs>

//             <div className="space-y-2">
//               <div className="flex items-center gap-2">
//                 <Label htmlFor="method-select" className="text-sm font-medium">
//                   Select method
//                 </Label>
//                 <div className="group relative">
//                   <Info className="w-4 h-4 text-gray-400 cursor-help" />
//                   <div className="absolute hidden group-hover:block bg-gray-500 text-white text-xs rounded p-2 w-80 -top-2 left-6 z-10">
//   <ul className="list-disc list-inside space-y-1">
//     <li>The Girolami method is simple and quick. It can be used for any compounds.</li>
//     <li>The Immirzi and Perini method is slightly more accurate and inconsiderably slower. It can be used for compounds with elements H, C, O, N, S, F, Cl, Br, I, Na, K and Rb. It is applicable only for hydrates, not for other solvates.</li>
//     <li>Refer to "International Journal of Pharmaceutics 355 (2008) 231–237" and the relevant literature cited within it for more information.</li>
//   </ul>
// </div>


//                 </div>
//               </div>
//               <Select value={selectedMethod} onValueChange={setSelectedMethod}>
//                 <SelectTrigger className="w-full h-10">
//                   <SelectValue placeholder="Choose prediction method" />
//                 </SelectTrigger>
//                 <SelectContent>
//                   <SelectItem value="girolami">Girolami</SelectItem>
//                   <SelectItem value="immirzi">Immirzi-Perini</SelectItem>
//                 </SelectContent>
//               </Select>
//             </div>

//             <Button
//               onClick={handleSubmit}
//               className="w-full bg-[var(--color-chart-2)] hover:bg-[var(--color-chart-3)] text-white"
//               disabled={isSubmitDisabled}
//             >
//               {isLoading ? <Loader2 className="w-4 h-4 mr-2 animate-spin" /> : <Play className="w-4 h-4 mr-2" />}
//               Submit
//             </Button>
//           </CardContent>
//         </Card>

//         {/* Results Section */}
//         <Card className="bg-white shadow-sm flex-1">
//           <CardHeader className="pb-4">
//             <CardTitle className="text-lg font-semibold ">
//               {selectedMethod && results.length > 0
//                 ? `${
//                     selectedMethod === "immirzi"
//                       ? "Immirzi-Perini"
//                       : selectedMethod.charAt(0).toUpperCase() + selectedMethod.slice(1)
//                   } Method: Results`
//                 : "Results"}
//             </CardTitle>
//           </CardHeader>
//           <CardContent className="flex items-center justify-center h-[calc(100%-4rem)]">
//             <ResultsTable results={results} isLoading={isLoading} inputType={activeTab === "single" ? "smiles" : "file"} />
//           </CardContent>
//         </Card>
//       </div>
//     </div>
//   );
// }



"use client";

import type React from "react";
import { useState, useEffect } from "react";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Textarea } from "@/components/ui/textarea";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Upload, Info, Play, Loader2, FileText } from "lucide-react";
import { ResultsTable, type PredictionResult } from "./ResultsTable";

export default function InputBar() {
  const [activeTab, setActiveTab] = useState("single");
  const [smilesInput, setSmilesInput] = useState("");
  const [multipleSmilesInput, setMultipleSmilesInput] = useState("");
  const [selectedMethod, setSelectedMethod] = useState("");
  const [uploadedFile, setUploadedFile] = useState<File | null>(null);
  const [results, setResults] = useState<PredictionResult[]>([]);
  const [isLoading, setIsLoading] = useState(false);

  // Sample data for the textarea
  const sampleSmilesData = `Phenacetin, CCOc1ccc(NC(C)=O)cc1
1, COc1ccc2c(c1)c(c(c(=O)n2C)CC(=O)O)Cl
A, CC(C(=O)O)c1ccc(cc1)C(=O)C`;

  // Sample data for single SMILES
  const singleSmilesSample = "CCOc1ccc(NC(C)=O)cc1";

  // Effect to clear results when switching tabs
  useEffect(() => {
    setResults([]); // Clear results whenever activeTab changes
  }, [activeTab]);

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
        const smilesArray = multipleSmilesInput.trim().split('\n');
        let csvContent = "compound,smiles\n";
        smilesArray.forEach(line => {
          const [compound, smiles] = line.split(',').map(s => s.trim());
          if (smiles) {
            csvContent += `${compound || `Molecule_${smilesArray.indexOf(line) + 1}`},${smiles}\n`;
          }
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
      setResults([{ smiles: "", compound: errorMessage }]); // Show error in compound field
    } finally {
      setIsLoading(false);
    }
  };

  const isSubmitDisabled = !selectedMethod || isLoading || (activeTab === "single" ? !smilesInput : (!uploadedFile && !multipleSmilesInput));

  const handleTextareaFocus = () => {
    setUploadedFile(null); // Clear uploaded file when focusing on textarea
  };

  const handleFileInputChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    setMultipleSmilesInput(""); // Clear textarea when uploading a file
    setUploadedFile(e.target.files?.[0] ?? null);
  };

  // Handle sample data loading for textarea
  const handleSampleDataLoad = () => {
    setMultipleSmilesInput(sampleSmilesData);
    setUploadedFile(null); // Clear uploaded file when using sample data
  };

  // Handle sample CSV loading
  const handleSampleCsvLoad = async () => {
    try {
      const response = await fetch('/molecules.csv');
      if (response.ok) {
        const csvContent = await response.text();
        const file = new File([csvContent], 'molecules.csv', { type: 'text/csv' });
        setUploadedFile(file);
        setMultipleSmilesInput(""); // Clear textarea when loading sample
      } else {
        alert('Sample file not found. Please ensure molecules.csv is in the public directory.');
      }
    } catch (error) {
      alert('Error loading sample file.');
    }
  };

  // Handle sample data loading for single SMILES
  const handleSingleSmilesSample = () => {
    setSmilesInput(singleSmilesSample);
  };

  return (
    <div className="max-w-screen-2xl mx-auto p-6">
      <div className="flex gap-6 h-[calc(100vh-8rem)]">
        {/* Input Section */}
        <Card className="bg-white shadow-sm w-[28%] flex-shrink-0">
          <CardHeader className="pb-1">
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
                <div className="flex items-center gap-2">
                  <Label htmlFor="smiles-input" className="text-sm font-medium">
                    SMILES
                  </Label>
                  <Button
                    variant="outline"
                    size="sm"
                    onClick={handleSingleSmilesSample}
                    className="px-2 ml-40 py-1 h-6 text-xs text-blue-600 border-blue-200 hover:bg-blue-50"
                  >
                    <FileText className="w-3 h-3 mr-1" />
                    Try Sample
                  </Button>
                </div>
                <Textarea
                  id="smiles-input"
                  placeholder="Enter single SMILES string"
                  value={smilesInput}
                  onChange={(e) => setSmilesInput(e.target.value)}
                  className="min-h-[100px] resize-none"
                />
              </TabsContent>
              <TabsContent value="multiple" className="mt-4 space-y-4">
                <div className="space-y-2">
                  <div className="flex items-center gap-2">
                    <Label className="text-sm font-medium">Upload csv</Label>
                    <div className="group relative">
                      <Info className="w-4 h-4 text-gray-400 " />
                      <span className="absolute hidden group-hover:block bg-gray-500 text-white text-xs rounded p-2 w-64 -top-2 left-6 z-10">
                        The compound IDs could be their generic names or any other custom identifiers given by the user.
                      </span>
                    </div>
                    <Button
                      variant="outline"
                      size="sm"
                      onClick={handleSampleCsvLoad}
                      className="px-2 ml-26 py-1 h-6 text-xs text-blue-600 border-blue-200 hover:bg-blue-50"
                    >
                      <FileText className="w-3 h-3 mr-1" />
                      Try Sample
                    </Button>
                  </div>
                  <div className="border-2 border-dashed border-gray-300 rounded-lg p-6 text-center bg-gray-50 relative">
                    <Upload className="w-10 h-5 text-gray-400 mx-auto mb-2 " />
                    <p className="text-xs font-medium text-gray-700">Click to upload or drag CSV file with multiple SMILES along with their compound IDs (new line separated).

</p>
                    <Input
                      type="file"
                      accept=".csv"
                      onChange={handleFileInputChange}
                      className="hidden"
                      id="file-upload"
                    />
                    <Label htmlFor="file-upload" className="absolute inset-0 cursor-pointer" />
                  </div>
                  {uploadedFile && <p className="text-sm text-green-600">File: {uploadedFile.name}</p>}
                </div>
                <div className="space-y-2">
                  <div className="flex items-center">
                    <Label htmlFor="multiple-smiles-input" className="text-sm font-medium">
                      Or enter multiple SMILES (newline separated)
                    </Label>
                    <div className="group relative ml-2">
                      <Info className="w-4 h-4 text-gray-400 " />
                      <span className="absolute hidden group-hover:block bg-gray-500 text-white text-xs rounded p-2 w-64 -top-2 left-6 z-10">
                        The compound IDs could be their generic names or any other custom identifiers given by the user.
                      </span>
                    </div>
                  </div>
                  <div className="flex items-center justify-between mt-2">
                    <Button
                      variant="outline"
                      size="sm"
                      onClick={handleSampleDataLoad}
                      className="px-2 py-1 h-6 text-xs text-blue-600 border-blue-200 hover:bg-blue-50"
                    >
                      <FileText className="w-3 h-3 mr-1" />
                      Try Sample
                    </Button>
                  </div>
                  <Textarea
                    id="multiple-smiles-input"
                    placeholder={`Phenacetin, CCOc1ccc(NC(C)=O)cc1
AGAN_999, Cc1cc(ccc1OCCCCCCN(S(=O)(=O)c1ccc(cc1)C(C)(C)C)C)C)S(=O)(=O)N(C)C
1, COc1ccc2c(c1)c(c(c(=O)n2C)CC(=O)O)Cl
A, CC(C(=O)O)c1ccc(cc1)C(=O)C`}
                    value={multipleSmilesInput}
                    onChange={(e) => setMultipleSmilesInput(e.target.value)}
                    onFocus={handleTextareaFocus}
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
                <div className="group relative">
                  <Info className="w-4 h-4 text-gray-400 " />
                  <div className="absolute hidden group-hover:block bg-gray-500 text-white text-xs rounded p-2 w-80 -top-2 left-6 z-10">
  <ul className="list-disc list-inside space-y-1">
    <li>The Girolami method is simple and quick. It can be used for any compound.</li>
    <li>The Immirzi-Perini method is slightly more accurate and marginally slower than the Girolami method. However, the Immirzi-Perini method can only be used for compounds containing the elements H, C, O, N, S, F, Cl, Br, I, Na, K, and Rb. Besides, it is applicable only to hydrates, not other solvates.</li>
    <li>Refer to the "International Journal of Pharmaceutics 355 (2008) 231–237" article and the relevant literature cited within it for more information on both methods.</li>
  </ul>
</div>
                </div>
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
            <ResultsTable results={results} isLoading={isLoading} inputType={activeTab === "single" ? "smiles" : "file"} />
          </CardContent>
        </Card>
      </div>
    </div>
  );
}