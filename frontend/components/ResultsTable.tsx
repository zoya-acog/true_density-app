"use client";

import type React from "react";
import { useState, useMemo, useCallback, useEffect } from "react";
import { Loader2, Info, Settings2, Filter, Download, X } from "lucide-react";
import { Button } from "@/components/ui/button";
import { Popover, PopoverContent, PopoverTrigger } from "@/components/ui/popover";
import { Checkbox } from "@/components/ui/checkbox";
import { Label } from "@/components/ui/label";
import Image from "next/image";

// Fallback to regular table if AG Grid is not available
let AgGridReact: any = null;
let ModuleRegistry: any = null;
let AllCommunityModule: any = null;
let agGridAvailable = false;

try {
  const agGrid = require("ag-grid-react");
  const agGridCommunity = require("ag-grid-community");
  
  AgGridReact = agGrid.AgGridReact;
  ModuleRegistry = agGridCommunity.ModuleRegistry;
  AllCommunityModule = agGridCommunity.AllCommunityModule;
  
  // Register AG Grid modules
  ModuleRegistry.registerModules([AllCommunityModule]);
  
  require("ag-grid-community/styles/ag-grid.css");
  require("ag-grid-community/styles/ag-theme-alpine.css");
  agGridAvailable = true;
} catch (error) {
  console.warn("AG Grid not available, falling back to regular table");
}

// Define types for AG Grid if available
type ColDef = any;
type GridReadyEvent = any;
type GridApi = any;
export interface PredictionResult {
  compound?: string; // Optional for single SMILES
  smiles: string;
  image?: string; // base64 encoded image
  density?: number;
  error_density?: number;
}

interface ResultsTableProps {
  results: PredictionResult[];
  isLoading: boolean;
  inputType?: string; // To determine single SMILES vs multiple
}

// Modal component for enlarged structure view
const StructureModal = ({ 
  isOpen, 
  onClose, 
  imageData, 
  smiles, 
  compound 
}: { 
  isOpen: boolean; 
  onClose: () => void; 
  imageData?: string; 
  smiles?: string; 
  compound?: string; 
}) => {
  if (!isOpen || !imageData) return null;

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center">
      {/* Backdrop */}
      <div 
        className="absolute inset-0 bg-black bg-opacity-50 backdrop-blur-sm"
        onClick={onClose}
      />
      
      {/* Modal content */}
      <div className="relative bg-white rounded-lg shadow-2xl max-w-2xl max-h-[90vh] overflow-hidden z-10">
        {/* Header */}
        <div className="flex items-center justify-between p-6 border-b border-gray-200">
          <div>
            <h2 className="text-xl font-semibold text-gray-900">Molecular Structure</h2>
            {compound && (
              <p className="text-sm text-gray-600 mt-1">{compound}</p>
            )}
          </div>
          <Button
            variant="ghost"
            size="sm"
            onClick={onClose}
            className="h-8 w-8 p-0 hover:bg-gray-100"
          >
            <X className="h-4 w-4" />
          </Button>
        </div>

        {/* Image content */}
        <div className="p-6 flex flex-col items-center">
          <div className="bg-white border border-gray-200 rounded-lg p-4 mb-4">
            <Image
              src={`data:image/png;base64,${imageData}`}
              alt={`Structure of ${smiles || 'compound'}`}
              width={400}
              height={400}
              className="object-contain"
              priority
            />
          </div>
          
          {/* SMILES string */}
          {smiles && (
            <div className="w-full">
              <Label className="text-sm font-medium text-gray-700">SMILES</Label>
              <div className="mt-1 p-3 bg-gray-50 rounded-md border">
                <code className="text-sm font-mono text-gray-800 break-all">
                  {smiles}
                </code>
              </div>
            </div>
          )}
        </div>
      </div>
    </div>
  );
};

// Fallback table component if AG Grid is not available
const FallbackTable = ({ 
  results, 
  visibleColumns, 
  onImageClick,
  inputType
}: { 
  results: PredictionResult[], 
  visibleColumns: any,
  onImageClick: (imageData: string, smiles: string, compound: string) => void,
  inputType?: string
}) => {
  const [currentPage, setCurrentPage] = useState(1);
  const [pageSize, setPageSize] = useState(10);
  const [sortField, setSortField] = useState<keyof PredictionResult | null>(null);
  const [sortDirection, setSortDirection] = useState<'asc' | 'desc'>('asc');
  const [filters, setFilters] = useState({
    compound: '',
    smiles: '',
    density: '',
    error_density: '',
  });

  // Filter and sort data
  const filteredAndSortedData = useMemo(() => {
    let filtered = results.filter(row => {
      return (
        (!filters.compound || row.compound?.toLowerCase().includes(filters.compound.toLowerCase())) &&
        (!filters.smiles || row.smiles?.toLowerCase().includes(filters.smiles.toLowerCase())) &&
        (!filters.density || (row.density?.toString().includes(filters.density))) &&
        (!filters.error_density || (row.error_density?.toString().includes(filters.error_density)))
      );
    });

    if (sortField) {
      filtered.sort((a, b) => {
        const aVal = a[sortField];
        const bVal = b[sortField];
        
        if (aVal === undefined && bVal === undefined) return 0;
        if (aVal === undefined) return 1;
        if (bVal === undefined) return -1;
        
        const comparison = aVal < bVal ? -1 : aVal > bVal ? 1 : 0;
        return sortDirection === 'desc' ? -comparison : comparison;
      });
    }

    return filtered;
  }, [results, filters, sortField, sortDirection]);

  // Paginate data
  const paginatedData = useMemo(() => {
    const startIndex = (currentPage - 1) * pageSize;
    return filteredAndSortedData.slice(startIndex, startIndex + pageSize);
  }, [filteredAndSortedData, currentPage, pageSize]);

  const totalPages = Math.ceil(filteredAndSortedData.length / pageSize);

  const handleSort = (field: keyof PredictionResult) => {
    if (sortField === field) {
      setSortDirection(sortDirection === 'asc' ? 'desc' : 'asc');
    } else {
      setSortField(field);
      setSortDirection('asc');
    }
  };

  return (
    <div className="space-y-4 mt-0">
      {/* Filters */}
      <div className="grid grid-cols-3 gap-4 p-4 bg-gray-50 rounded-lg">
        {inputType !== 'smiles' && visibleColumns.compound && (
          <div>
            <Label className="text-xs font-medium text-gray-500">Filter Compound</Label>
            <input
              type="text"
              className="w-full mt-1 px-3 py-2 text-sm border border-gray-300 rounded-md"
              placeholder="Search Compound..."
              value={filters.compound}
              onChange={(e) => setFilters(prev => ({ ...prev, compound: e.target.value }))}
            />
          </div>
        )}
        <div>
          <Label className="text-xs font-medium text-gray-500">Filter SMILES</Label>
          <input
            type="text"
            className="w-full mt-1 px-3 py-2 text-sm border border-gray-300 rounded-md"
            placeholder="Search SMILES..."
            value={filters.smiles}
            onChange={(e) => setFilters(prev => ({ ...prev, smiles: e.target.value }))}
          />
        </div>
      </div>

      {/* Table */}
      <div className="border border-gray-200 rounded-lg overflow-hidden">
        <table className="w-full">
          <thead className="bg-gray-50">
            <tr>
              {inputType !== 'smiles' && visibleColumns.compound && (
                <th
                  className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider cursor-pointer hover:bg-gray-100"
                  onClick={() => handleSort('compound')}
                >
                  Compound {sortField === 'compound' && (sortDirection === 'asc' ? '↑' : '↓')}
                </th>
              )}
              {visibleColumns.smiles && (
                <th 
                  className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider cursor-pointer hover:bg-gray-100"
                  onClick={() => handleSort('smiles')}
                >
                  SMILES {sortField === 'smiles' && (sortDirection === 'asc' ? '↑' : '↓')}
                </th>
              )}
              {visibleColumns.density && (
                <th 
                  className="px-6 py-3 text-right text-xs font-medium text-gray-500 uppercase tracking-wider cursor-pointer hover:bg-gray-100"
                  onClick={() => handleSort('density')}
                >
                  True Density(g/cm³) {sortField === 'density' && (sortDirection === 'asc' ? '↑' : '↓')}
                </th>
              )}
              {visibleColumns.error_density && (
                <th 
                  className="px-6 py-3 text-right text-xs font-medium text-gray-500 uppercase tracking-wider cursor-pointer hover:bg-gray-100"
                  onClick={() => handleSort('error_density')}
                >
                  <div className="text-center">
                    <div>True Density</div>
                    <div>with error correction</div>
                  </div>
                  {sortField === 'error_density' && (sortDirection === 'asc' ? '↑' : '↓')}
                </th>
              )}
              {visibleColumns.image && (
                <th className="px-6 py-3 text-center text-xs font-medium text-gray-500 uppercase tracking-wider">
                  2D - Structure
                </th>
              )}
            </tr>
          </thead>
          <tbody className="bg-white divide-y divide-gray-200">
            {paginatedData.map((result, index) => (
              <tr key={index} className="hover:bg-gray-50">
                {inputType !== 'smiles' && visibleColumns.compound && (
                  <td className="px-6 py-4 whitespace-nowrap text-sm font-mono text-gray-900">
                    <div className="max-w-xs truncate" title={result.compound}>
                      {result.compound}
                    </div>
                  </td>
                )}
                {visibleColumns.smiles && (
                  <td className="px-6 py-4 whitespace-nowrap text-sm font-mono text-gray-900">
                    <div className="max-w-xs truncate" title={result.smiles}>
                      {result.smiles}
                    </div>
                  </td>
                )}
                {visibleColumns.density && (
                  <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-900 text-right font-mono">
                    {result.density?.toFixed(3) ?? '-'}
                  </td>
                )}
                {visibleColumns.error_density && (
                  <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-900 text-right font-mono">
                    {result.error_density?.toFixed(3) ?? '-'}
                  </td>
                )}
                {visibleColumns.image && (
                  <td className="px-6 py-4 whitespace-nowrap text-center">
                    {result.image ? (
                      <div className="flex justify-center">
                        <Image
                          src={`data:image/png;base64,${result.image}`}
                          alt={`Structure of ${result.smiles}`}
                          width={60}
                          height={60}
                          className="object-contain cursor-pointer hover:opacity-80 transition-opacity"
                          onClick={() => onImageClick(result.image!, result.smiles, result.compound || '')}
                        />
                      </div>
                    ) : (
                      <div className="text-gray-400 text-xs">No image</div>
                    )}
                  </td>
                )}
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      {/* Pagination */}
      <div className="flex items-center justify-between">
        <div className="flex items-center gap-2">
          <Label className="text-sm text-gray-600">Rows per page:</Label>
          <select
            className="border border-gray-300 rounded px-2 py-1 text-sm"
            value={pageSize}
            onChange={(e) => {
              setPageSize(Number(e.target.value));
              setCurrentPage(1);
            }}
          >
            <option value={10}>10</option>
            <option value={25}>25</option>
            <option value={50}>50</option>
            <option value={100}>100</option>
          </select>
        </div>
        
        <div className="flex items-center gap-2">
          <span className="text-sm text-gray-600">
            {((currentPage - 1) * pageSize) + 1}-{Math.min(currentPage * pageSize, filteredAndSortedData.length)} of {filteredAndSortedData.length}
          </span>
          <div className="flex gap-1">
            <Button
              variant="outline"
              size="sm"
              disabled={currentPage === 1}
              onClick={() => setCurrentPage(prev => prev - 1)}
            >
              Previous
            </Button>
            <Button
              variant="outline"
              size="sm"
              disabled={currentPage === totalPages}
              onClick={() => setCurrentPage(prev => prev + 1)}
            >
              Next
            </Button>
          </div>
        </div>
      </div>
    </div>
  );
};

// Custom header component for multiline headers with sorting support
const MultilineHeader = (props: any) => {
  const { displayName, enableSorting, progressSort, column } = props;

  const onHeaderClick = (event: React.MouseEvent) => {
    if (enableSorting) {
      progressSort(event.shiftKey);
    }
  };

  const sort = column.getSort();
  const sortIcon = sort ? (sort === 'asc' ? '↑' : '↓') : null;

  const lines = displayName.split('|');

  return (
    <div
      style={{ 
        textAlign: 'center', 
        lineHeight: '1.2', 
        cursor: enableSorting ? 'pointer' : 'default',
        display: 'flex',
        flexDirection: 'column',
        justifyContent: 'center',
        height: '100%',
        width: '100%'
      }}
      onClick={onHeaderClick}
    >
      {lines.map((line: string, index: React.Key | null | undefined) => (
        <div key={index}>{line.trim()}</div>
      ))}
      {sortIcon && <span style={{ fontSize: '12px', marginTop: '2px' }}>{sortIcon}</span>}
    </div>
  );
};

// Custom cell renderer for structure images with click handler
const ImageCellRenderer = ({ 
  value, 
  data, 
  onImageClick 
}: { 
  value?: string; 
  data: PredictionResult;
  onImageClick: (imageData: string, smiles: string, compound: string) => void;
}) => {
  if (!value) return <div className="text-gray-400 text-sm">No image</div>;
  
  return (
    <div className="flex justify-center items-center h-full">
      <Image
        src={`data:image/png;base64,${value}`}
        alt="Molecular structure"
        width={60}
        height={60}
        className="object-contain cursor-pointer hover:opacity-80 transition-opacity hover:shadow-lg rounded"
        onClick={() => onImageClick(value, data.smiles, data.compound || '')}
      />
    </div>
  );
};

// Custom cell renderer for numeric values with formatting
const NumericCellRenderer = ({ value }: { value?: number }) => {
  if (value === undefined || value === null) {
    return <div className="text-gray-400">-</div>;
  }
  return <div className="text-right font-mono">{value.toFixed(3)}</div>;
};

export function ResultsTable({ results, isLoading, inputType }: ResultsTableProps) {
  const [gridApi, setGridApi] = useState<GridApi | null>(null);
  const isSingleSmiles = inputType === 'smiles';
  const [visibleColumns, setVisibleColumns] = useState({
    compound: !isSingleSmiles, // Ensure compound is visible by default for file inputs
    smiles: true,
    density: true,
    error_density: true,
    image: true,
  });
  
  // Modal state
  const [modalState, setModalState] = useState({
    isOpen: false,
    imageData: '',
    smiles: '',
    compound: ''
  });

  // Handle image click
  const handleImageClick = useCallback((imageData: string, smiles: string, compound: string) => {
    setModalState({
      isOpen: true,
      imageData,
      smiles,
      compound
    });
  }, []);

  // Close modal
  const closeModal = useCallback(() => {
    setModalState(prev => ({
      ...prev,
      isOpen: false
    }));
  }, []);

  // Close modal on Escape key
  useEffect(() => {
    const handleEscape = (event: KeyboardEvent) => {
      if (event.key === 'Escape') {
        closeModal();
      }
    };

    if (modalState.isOpen) {
      document.addEventListener('keydown', handleEscape);
      // Prevent body scroll when modal is open
      document.body.style.overflow = 'hidden';
    } else {
      document.body.style.overflow = 'unset';
    }

    return () => {
      document.removeEventListener('keydown', handleEscape);
      document.body.style.overflow = 'unset';
    };
  }, [modalState.isOpen, closeModal]);

  // Column definitions
  const columnDefs: ColDef[] = useMemo(() => {
    const columns = [
      {
        headerName: "SMILES",
        sortable: false,
        field: "smiles",
        flex: 1,
        minWidth: 200,
        filter: "agTextColumnFilter",
        floatingFilter: true,
        cellClass: "font-mono text-sm",
        tooltipField: "smiles",
        hide: !visibleColumns.smiles,
        
      },
      {
        headerName: "True Density (g/cm³)",
        field: "density",
        flex: 1,
        minWidth: 150,
        filter: "agNumberColumnFilter",
        floatingFilter: true,
        cellRenderer: NumericCellRenderer,
        hide: !visibleColumns.density,
        comparator: (valueA: number, valueB: number) => {
          if (valueA == null && valueB == null) return 0;
          if (valueA == null) return -1;
          if (valueB == null) return 1;
          return valueA - valueB;
        },
        cellStyle: { display: 'flex', alignItems: 'center', justifyContent: 'center' },
      },
      {
        headerName: "True Density|with error correction|(g/cm³)",
        field: "error_density",
        flex: 1,
        minWidth: 170,
        filter: "agNumberColumnFilter",
        floatingFilter: true,
        cellRenderer: NumericCellRenderer,
        hide: !visibleColumns.error_density,
        headerComponent: MultilineHeader,
        comparator: (valueA: number, valueB: number) => {
          if (valueA == null && valueB == null) return 0;
          if (valueA == null) return -1;
          if (valueB == null) return 1;
          return valueA - valueB;
        },
        cellStyle: { display: 'flex', alignItems: 'center', justifyContent: 'center' },
        
      },
      {
        headerName: "2D - Structure",
        field: "image",
        flex: 1,
        minWidth: 120,
        sortable: false,
        filter: false,
        cellRenderer: (params: any) => ImageCellRenderer({
          value: params.value,
          data: params.data,
          onImageClick: handleImageClick
        }),
        cellStyle: { display: 'flex', alignItems: 'center', justifyContent: 'center' },
        hide: !visibleColumns.image
      }
    ];

    if (!isSingleSmiles) {
      columns.unshift({
        headerName: "Compound",
        field: "compound",
        sortable: true, // Added to fix the TypeScript error
        flex: 1,
        minWidth: 150,
        filter: "agTextColumnFilter",
        floatingFilter: true,
        cellClass: "font-mono text-sm",
        tooltipField: "compound",
        hide: false // Ensure compound column is visible by default for file inputs
      });
    }

    return columns;
  }, [visibleColumns, handleImageClick, isSingleSmiles]);

  // Default column configuration
  const defaultColDef: ColDef = useMemo(() => ({
    sortable: true,
    resizable: true,
    suppressSizeToFit: false,
  }), []);

  // Grid ready handler
  const onGridReady = useCallback((params: GridReadyEvent) => {
    setGridApi(params.api);
  }, []);

  // Toggle column visibility
  const toggleColumn = useCallback((columnKey: keyof typeof visibleColumns) => {
    setVisibleColumns(prev => ({
      ...prev,
      [columnKey]: !prev[columnKey]
    }));
  }, []);

  // Clear all filters
  const clearFilters = useCallback(() => {
    if (agGridAvailable && gridApi) {
      gridApi.setFilterModel(null);
    }
  }, [gridApi]);

  // Export data as CSV
  const exportToCsv = useCallback(() => {
    if (agGridAvailable && gridApi) {
      gridApi.exportDataAsCsv({
        fileName: 'prediction_results.csv'
      });
    } else {
      // Fallback CSV export
      const headers = Object.keys(visibleColumns).filter(key => visibleColumns[key as keyof typeof visibleColumns] && (key !== 'compound' || !isSingleSmiles));
      const csvContent = [
        headers.join(','),
        ...results.map(row => 
          headers.map(header => {
            const value = row[header as keyof PredictionResult];
            return typeof value === 'string' ? `"${value}"` : (value ?? '');
          }).join(',')
        )
      ].join('\n');
      
      const blob = new Blob([csvContent], { type: 'text/csv' });
      const url = window.URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = 'prediction_results.csv';
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      window.URL.revokeObjectURL(url);
    }
  }, [gridApi, results, visibleColumns, isSingleSmiles]);

  if (isLoading) {
    return (
      <div className="flex flex-col items-center justify-center h-full text-gray-500">
        <Loader2 className="w-12 h-12 animate-spin text-gray-400" />
        <p className="mt-4">Processing, please wait...</p>
      </div>
    );
  }

  if (results.length === 0) {
    return (
      <div className="text-center text-gray-500 space-y-4">
        <div className="w-16 h-16 bg-gray-200 rounded-full flex items-center justify-center mx-auto">
          <Info className="w-8 h-8 text-gray-400" />
        </div>
        <p className="text-sm">
          Upload a CSV file or enter SMILES and select a model to see results.
        </p>
      </div>
    );
  }

  return (
    <>
      <div className="w-full h-full flex flex-col">
        <div className="flex items-center justify-between p-4 border-b border-gray-200">
          <div className="flex items-center gap-2">
          </div>
          <div className="flex items-center gap-2">
            <Button
              variant="outline"
              size="sm"
              onClick={clearFilters}
              className="flex items-center gap-2"
            >
              <Filter className="w-4 h-4" />
              Clear Filters
            </Button>
            
            <Popover>
              <PopoverTrigger asChild>
                <Button variant="outline" size="sm" className="flex items-center gap-2">
                  <Settings2 className="w-4 h-4" />
                  Columns
                </Button>
              </PopoverTrigger>
              <PopoverContent className="w-56" align="end">
                <div className="space-y-4">
                  <h4 className="font-medium text-sm">Toggle Columns</h4>
                  <div className="space-y-3">
                    {Object.entries(visibleColumns).map(([key, visible]) => (
                      (key !== 'compound' || !isSingleSmiles) && (
                      <div key={key} className="flex items-center space-x-2">
                        <Checkbox
                          id={key}
                          checked={visible}
                          onCheckedChange={() => toggleColumn(key as keyof typeof visibleColumns)}
                        />
                        <Label htmlFor={key} className="text-sm font-normal capitalize">
                          {key.replace('_', ' ')}
                        </Label>
                      </div>
                    )))}
                  </div>
                </div>
              </PopoverContent>
            </Popover>

            <Button
              variant="outline"
              size="sm"
              onClick={exportToCsv}
              className="flex items-center gap-2"
            >
              <Download className="w-4 h-4" />
              Export CSV
            </Button>
          </div>
        </div>

        <div className="flex-1 min-h-0">
          <style jsx>{`
            .ag-theme-alpine {
              --ag-background-color: white;
              --ag-header-background-color: #f9fafb;
              --ag-header-foreground-color: #6b7280;
              --ag-border-color: #b8b8bcff;
              --ag-row-border-color: #e5e7eb;
              --ag-header-height: 60px;
              --ag-row-height: 78px;
              --ag-font-size: 14px;
              --ag-font-family: ui-sans-serif, system-ui, sans-serif, "Apple Color Emoji", "Segoe UI Emoji", "Segoe UI Symbol", "Noto Color Emoji";
            }
            .ag-theme-alpine .ag-header {
              border-bottom: 1px solid #e5e7eb;
              border-top: none;
            }
            .ag-theme-alpine .ag-header-cell {
              padding: 12px 24px;
              font-size: 11px;
              font-weight: 500;
              color: #6b7280;
              text-transform: uppercase;
              letter-spacing: 0.05em;
              border-right: none;
              background: #f9fafb;
              font-family: ui-sans-serif, system-ui, sans-serif, "Apple Color Emoji", "Segoe UI Emoji", "Segoe UI Symbol", "Noto Color Emoji";
              display: flex;
              alignItems: center;
              justifyContent: center;
            }
            .ag-theme-alpine .ag-header-cell-text {
              overflow: visible;
              white-space: pre-line;
              text-align: center;
            }
            .ag-theme-alpine .multiline-header .ag-header-cell-text {
              line-height: 1.2;
              white-space: pre-line;
              text-align: center;
            }
            .ag-theme-alpine .ag-cell {
              padding: 20px 24px;
              border-right: none;
              border-bottom: none;
              display: flex;
              alignItems: center;
              background: white;
              font-size: 14px;
              font-family: ui-sans-serif, system-ui, sans-serif, "Apple Color Emoji", "Segoe UI Emoji", "Segoe UI Symbol", "Noto Color Emoji";
            }
            .ag-theme-alpine .ag-row {
              border-bottom: 1px solid #e5e7eb;
              background: white;
            }
            .ag-theme-alpine .ag-row:hover {
              background-color: #f9fafb;
            }
            .ag-theme-alpine .ag-root-wrapper {
              border: none;
              background: white;
            }
            .ag-theme-alpine .ag-root {
              border: none;
            }
            .ag-theme-alpine .ag-body-viewport {
              background: white;
            }
            .ag-theme-alpine .ag-row:last-child {
              border-bottom: none;
            }
            .ag-theme-alpine .ag-body-horizontal-scroll-viewport {
              overflow: visible;
            }
            .ag-theme-alpine .ag-floating-filter-input {
              height: 32px;
              font-size: 12px;
            }
            .ag-theme-alpine .ag-floating-filter-wrapper {
              padding: 8px 12px;
            }
          `}</style>
          
          <div className="ag-theme-alpine h-full">
            <AgGridReact
              rowData={results}
              columnDefs={columnDefs}
              defaultColDef={defaultColDef}
              onGridReady={onGridReady}
              domLayout="normal"
              headerHeight={60}
              rowHeight={88}
              suppressRowHoverHighlight={false}
              pagination={true}
              paginationPageSize={10}
              paginationPageSizeSelector={[10, 25, 50, 100]}
              animateRows={true}
              enableCellTextSelection={true}
              suppressCopyRowsToClipboard={false}
              suppressMovableColumns={false}
              suppressColumnVirtualisation={true}
              suppressRowVirtualisation={false}
              theme="legacy"
            />
          </div>
        </div>
      </div>

      <StructureModal
        isOpen={modalState.isOpen}
        onClose={closeModal}
        imageData={modalState.imageData}
        smiles={modalState.smiles}
        compound={modalState.compound}
      />
    </>
  );
}