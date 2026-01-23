'use client';

import { useState, useMemo } from 'react';
import { cn } from '@/lib/utils';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Checkbox } from '@/components/ui/checkbox';
import { Badge } from '@/components/ui/badge';
import { Input } from '@/components/ui/input';
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from '@/components/ui/table';
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from '@/components/ui/select';
import {
  List,
  Download,
  CheckCircle,
  AlertTriangle,
  XCircle,
  ArrowUpDown,
  Filter,
  FileText,
} from 'lucide-react';
import type { PipelineDesign } from '@/lib/store';

interface PipelineResultsViewProps {
  results: PipelineDesign[];
  onExportFasta?: (selectedIds: string[], includeReview: boolean) => void;
  onSelectForAF3?: (designs: PipelineDesign[]) => void;
}

type SortField = 'name' | 'config' | 'plddt' | 'ptm' | 'pae' | 'status';
type SortDirection = 'asc' | 'desc';
type StatusFilter = 'all' | 'pass' | 'review' | 'fail';

export function PipelineResultsView({
  results,
  onExportFasta,
  onSelectForAF3,
}: PipelineResultsViewProps) {
  const [sortField, setSortField] = useState<SortField>('plddt');
  const [sortDirection, setSortDirection] = useState<SortDirection>('desc');
  const [statusFilter, setStatusFilter] = useState<StatusFilter>('all');
  const [searchQuery, setSearchQuery] = useState('');
  const [selectedIds, setSelectedIds] = useState<Set<string>>(new Set());

  // Filter and sort results
  const filteredResults = useMemo(() => {
    let filtered = [...results];

    // Apply status filter
    if (statusFilter !== 'all') {
      filtered = filtered.filter((d) => d.status === statusFilter);
    }

    // Apply search filter
    if (searchQuery) {
      const query = searchQuery.toLowerCase();
      filtered = filtered.filter(
        (d) =>
          d.name.toLowerCase().includes(query) ||
          d.config.toLowerCase().includes(query)
      );
    }

    // Sort
    filtered.sort((a, b) => {
      let comparison = 0;

      switch (sortField) {
        case 'name':
          comparison = a.name.localeCompare(b.name);
          break;
        case 'config':
          comparison = a.config.localeCompare(b.config);
          break;
        case 'plddt':
          comparison = a.plddt - b.plddt;
          break;
        case 'ptm':
          comparison = a.ptm - b.ptm;
          break;
        case 'pae':
          comparison = a.pae - b.pae;
          break;
        case 'status':
          const statusOrder = { pass: 0, review: 1, fail: 2 };
          comparison = statusOrder[a.status] - statusOrder[b.status];
          break;
      }

      return sortDirection === 'asc' ? comparison : -comparison;
    });

    return filtered;
  }, [results, sortField, sortDirection, statusFilter, searchQuery]);

  // Toggle sort
  const toggleSort = (field: SortField) => {
    if (sortField === field) {
      setSortDirection((d) => (d === 'asc' ? 'desc' : 'asc'));
    } else {
      setSortField(field);
      setSortDirection('desc');
    }
  };

  // Toggle selection
  const toggleSelection = (id: string) => {
    const newSelection = new Set(selectedIds);
    if (newSelection.has(id)) {
      newSelection.delete(id);
    } else {
      newSelection.add(id);
    }
    setSelectedIds(newSelection);
  };

  // Select all visible
  const selectAllVisible = () => {
    const newSelection = new Set(selectedIds);
    filteredResults.forEach((d) => newSelection.add(d.name));
    setSelectedIds(newSelection);
  };

  // Clear selection
  const clearSelection = () => {
    setSelectedIds(new Set());
  };

  // Get status badge
  const getStatusBadge = (status: PipelineDesign['status']) => {
    switch (status) {
      case 'pass':
        return (
          <Badge variant="default" className="bg-green-500">
            <CheckCircle className="w-3 h-3 mr-1" />
            Pass
          </Badge>
        );
      case 'review':
        return (
          <Badge variant="secondary" className="bg-yellow-500/20 text-yellow-600">
            <AlertTriangle className="w-3 h-3 mr-1" />
            Review
          </Badge>
        );
      case 'fail':
        return (
          <Badge variant="destructive">
            <XCircle className="w-3 h-3 mr-1" />
            Fail
          </Badge>
        );
    }
  };

  // Count by status
  const statusCounts = useMemo(() => {
    const counts = { all: results.length, pass: 0, review: 0, fail: 0 };
    results.forEach((d) => {
      counts[d.status]++;
    });
    return counts;
  }, [results]);

  return (
    <Card>
      <CardHeader>
        <CardTitle className="text-sm flex items-center justify-between">
          <div className="flex items-center gap-2">
            <List className="w-4 h-4" />
            Pipeline Results
            <Badge variant="secondary">{results.length} designs</Badge>
          </div>
          <div className="flex gap-2">
            {selectedIds.size > 0 && (
              <Button
                variant="outline"
                size="sm"
                onClick={() =>
                  onExportFasta?.(Array.from(selectedIds), false)
                }
              >
                <Download className="w-4 h-4 mr-2" />
                Export Selected ({selectedIds.size})
              </Button>
            )}
            {onSelectForAF3 && selectedIds.size > 0 && (
              <Button
                variant="default"
                size="sm"
                onClick={() => {
                  const selected = results.filter((d) =>
                    selectedIds.has(d.name)
                  );
                  onSelectForAF3(selected);
                }}
              >
                <FileText className="w-4 h-4 mr-2" />
                Send to AF3
              </Button>
            )}
          </div>
        </CardTitle>
      </CardHeader>
      <CardContent className="space-y-4">
        {/* Filters */}
        <div className="flex gap-4 items-center">
          <Input
            placeholder="Search by name or config..."
            value={searchQuery}
            onChange={(e) => setSearchQuery(e.target.value)}
            className="max-w-xs"
          />

          <Select value={statusFilter} onValueChange={(v) => setStatusFilter(v as StatusFilter)}>
            <SelectTrigger className="w-[180px]">
              <div className="flex items-center gap-2">
                <Filter className="w-4 h-4" />
                <SelectValue placeholder="Filter by status" />
              </div>
            </SelectTrigger>
            <SelectContent>
              <SelectItem value="all">
                All ({statusCounts.all})
              </SelectItem>
              <SelectItem value="pass">
                <div className="flex items-center gap-2">
                  <CheckCircle className="w-4 h-4 text-green-500" />
                  Pass ({statusCounts.pass})
                </div>
              </SelectItem>
              <SelectItem value="review">
                <div className="flex items-center gap-2">
                  <AlertTriangle className="w-4 h-4 text-yellow-500" />
                  Review ({statusCounts.review})
                </div>
              </SelectItem>
              <SelectItem value="fail">
                <div className="flex items-center gap-2">
                  <XCircle className="w-4 h-4 text-red-500" />
                  Fail ({statusCounts.fail})
                </div>
              </SelectItem>
            </SelectContent>
          </Select>

          <div className="flex-1" />

          <Button
            variant="ghost"
            size="sm"
            onClick={selectAllVisible}
            disabled={filteredResults.length === 0}
          >
            Select All Visible
          </Button>
          <Button
            variant="ghost"
            size="sm"
            onClick={clearSelection}
            disabled={selectedIds.size === 0}
          >
            Clear Selection
          </Button>
        </div>

        {/* Results Table */}
        <div className="rounded-md border">
          <Table>
            <TableHeader>
              <TableRow>
                <TableHead className="w-10">
                  <Checkbox
                    checked={
                      filteredResults.length > 0 &&
                      filteredResults.every((d) => selectedIds.has(d.name))
                    }
                    onCheckedChange={(checked) => {
                      if (checked) {
                        selectAllVisible();
                      } else {
                        clearSelection();
                      }
                    }}
                  />
                </TableHead>
                <SortableHeader
                  field="name"
                  label="Name"
                  currentField={sortField}
                  direction={sortDirection}
                  onSort={toggleSort}
                />
                <SortableHeader
                  field="config"
                  label="Config"
                  currentField={sortField}
                  direction={sortDirection}
                  onSort={toggleSort}
                />
                <SortableHeader
                  field="plddt"
                  label="pLDDT"
                  currentField={sortField}
                  direction={sortDirection}
                  onSort={toggleSort}
                />
                <SortableHeader
                  field="ptm"
                  label="pTM"
                  currentField={sortField}
                  direction={sortDirection}
                  onSort={toggleSort}
                />
                <SortableHeader
                  field="pae"
                  label="PAE"
                  currentField={sortField}
                  direction={sortDirection}
                  onSort={toggleSort}
                />
                <SortableHeader
                  field="status"
                  label="Status"
                  currentField={sortField}
                  direction={sortDirection}
                  onSort={toggleSort}
                />
              </TableRow>
            </TableHeader>
            <TableBody>
              {filteredResults.length === 0 ? (
                <TableRow>
                  <TableCell colSpan={7} className="text-center text-muted-foreground py-8">
                    No designs match the current filters
                  </TableCell>
                </TableRow>
              ) : (
                filteredResults.map((design) => (
                  <TableRow
                    key={design.name}
                    className={cn(
                      selectedIds.has(design.name) && 'bg-muted/50'
                    )}
                  >
                    <TableCell>
                      <Checkbox
                        checked={selectedIds.has(design.name)}
                        onCheckedChange={() => toggleSelection(design.name)}
                      />
                    </TableCell>
                    <TableCell className="font-mono text-sm">
                      {design.name}
                    </TableCell>
                    <TableCell className="text-muted-foreground text-sm">
                      {design.config}
                    </TableCell>
                    <TableCell>
                      <MetricCell
                        value={design.plddt}
                        format="percent"
                        goodThreshold={0.80}
                        warnThreshold={0.70}
                      />
                    </TableCell>
                    <TableCell>
                      <MetricCell
                        value={design.ptm}
                        format="percent"
                        goodThreshold={0.80}
                        warnThreshold={0.65}
                      />
                    </TableCell>
                    <TableCell>
                      <MetricCell
                        value={design.pae}
                        format="angstrom"
                        goodThreshold={5.0}
                        warnThreshold={10.0}
                        invertColors
                      />
                    </TableCell>
                    <TableCell>{getStatusBadge(design.status)}</TableCell>
                  </TableRow>
                ))
              )}
            </TableBody>
          </Table>
        </div>

        {/* Summary */}
        <div className="flex justify-between items-center text-sm text-muted-foreground">
          <span>
            Showing {filteredResults.length} of {results.length} designs
          </span>
          <span>{selectedIds.size} selected</span>
        </div>
      </CardContent>
    </Card>
  );
}

function SortableHeader({
  field,
  label,
  currentField,
  direction,
  onSort,
}: {
  field: SortField;
  label: string;
  currentField: SortField;
  direction: SortDirection;
  onSort: (field: SortField) => void;
}) {
  const isActive = field === currentField;

  return (
    <TableHead>
      <button
        className="flex items-center gap-1 hover:text-foreground transition-colors"
        onClick={() => onSort(field)}
      >
        {label}
        <ArrowUpDown
          className={cn(
            'w-3 h-3',
            isActive ? 'opacity-100' : 'opacity-50',
            isActive && direction === 'asc' && 'rotate-180'
          )}
        />
      </button>
    </TableHead>
  );
}

function MetricCell({
  value,
  format,
  goodThreshold,
  warnThreshold,
  invertColors = false,
}: {
  value: number;
  format: 'percent' | 'angstrom';
  goodThreshold: number;
  warnThreshold: number;
  invertColors?: boolean;
}) {
  let isGood: boolean;
  let isWarn: boolean;

  if (invertColors) {
    // For metrics where lower is better (PAE)
    isGood = value <= goodThreshold;
    isWarn = value <= warnThreshold && value > goodThreshold;
  } else {
    // For metrics where higher is better (pLDDT, pTM)
    isGood = value >= goodThreshold;
    isWarn = value >= warnThreshold && value < goodThreshold;
  }

  const displayValue =
    format === 'percent'
      ? `${(value * 100).toFixed(1)}%`
      : `${value.toFixed(1)} A`;

  return (
    <span
      className={cn(
        'font-medium',
        isGood && 'text-green-600 dark:text-green-400',
        isWarn && 'text-yellow-600 dark:text-yellow-400',
        !isGood && !isWarn && 'text-red-600 dark:text-red-400'
      )}
    >
      {displayValue}
    </span>
  );
}
