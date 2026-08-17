#!/bin/bash
#
# Example memory settings for a large-memory host (adjust as needed)
# Phase 2 annotation pipeline
#

# ============================================================================
# MEMORY ALLOCATION STRATEGY
# ============================================================================

# Example: 128 GB host (edit MAX_MEMORY_GB for your machine)
# Reserve 8GB for system processes
# Available for analysis: 120GB

# Memory allocation per tool:
MAX_MEMORY_GB=120
SYSTEM_RESERVE_GB=8

# Tool-specific memory limits
LAST_MEMORY_GB=30          # LAST alignment (most memory intensive)
SNPEFF_MEMORY_GB=80        # SnpEff annotation (Java heap)
PYTHON_MEMORY_GB=20        # Python scripts
BUFFER_MEMORY_GB=10        # Buffer for other processes

# ============================================================================
# JAVA OPTIONS FOR SNPEFF
# ============================================================================

# SnpEff Java memory settings
export JAVA_OPTS="-Xmx${SNPEFF_MEMORY_GB}g -Xms4g -XX:+UseG1GC -XX:MaxGCPauseMillis=200"

# ============================================================================
# OPTIONAL PROCESS LIMITS
# ============================================================================

ulimit -v $((MAX_MEMORY_GB * 1024 * 1024)) 2>/dev/null || true
ulimit -m $((MAX_MEMORY_GB * 1024)) 2>/dev/null || true

# ============================================================================
# OPTIONAL SYSTEM TUNING (disabled by default; requires root)
# ============================================================================

# echo 1 > /proc/sys/vm/drop_caches 2>/dev/null || true
# echo 1 > /proc/sys/vm/compact_memory 2>/dev/null || true

# ============================================================================
# MONITORING FUNCTIONS
# ============================================================================

check_memory_usage() {
    local used_memory=$(free -g | awk 'NR==2{print $3}')
    local available_memory=$(free -g | awk 'NR==2{print $7}')
    local memory_percent=$((used_memory * 100 / (used_memory + available_memory)))
    
    echo "Memory Status:"
    echo "  Used: ${used_memory}GB (${memory_percent}%)"
    echo "  Available: ${available_memory}GB"
    
    if [ "$used_memory" -gt "$MAX_MEMORY_GB" ]; then
        echo "  WARNING: Memory usage exceeds limit!"
        return 1
    fi
    
    return 0
}

monitor_memory() {
    local interval=${1:-60}  # Default 60 seconds
    
    while true; do
        if ! check_memory_usage; then
            echo "Memory limit exceeded. Consider stopping some processes."
            break
        fi
        sleep $interval
    done
}

# ============================================================================
# CHUNK PROCESSING FUNCTIONS
# ============================================================================

# Function to process large files in chunks
process_vcf_chunked() {
    local vcf_file="$1"
    local chunk_size="${2:-10000}"
    local output_dir="$3"
    
    # Create temporary directory for chunks
    local temp_dir="${output_dir}/chunks"
    mkdir -p "$temp_dir"
    
    # Split VCF into chunks
    local chunk_num=0
    local current_chunk=""
    local line_count=0
    
    while IFS= read -r line; do
        if [[ "$line" =~ ^# ]]; then
            # Header lines go to all chunks
            echo "$line" >> "${temp_dir}/chunk_${chunk_num}.vcf"
        else
            # Variant lines
            if [ $line_count -eq 0 ]; then
                # Start new chunk
                chunk_num=$((chunk_num + 1))
                current_chunk="${temp_dir}/chunk_${chunk_num}.vcf"
                # Copy header to new chunk
                head -n 1000 "$vcf_file" | grep "^#" > "$current_chunk"
            fi
            
            echo "$line" >> "$current_chunk"
            line_count=$((line_count + 1))
            
            if [ $line_count -ge $chunk_size ]; then
                line_count=0
            fi
        fi
    done < "$vcf_file"
    
    echo "$temp_dir"
}

# Clean up chunk files
cleanup_chunks() {
    local chunk_dir="$1"
    if [ -d "$chunk_dir" ]; then
        rm -rf "$chunk_dir"
    fi
}

# ============================================================================
# PARALLEL PROCESSING CONTROL
# ============================================================================

# Limit parallel processes based on memory
get_optimal_threads() {
    local tool="$1"
    local available_memory=$(free -g | awk 'NR==2{print $7}')
    
    case "$tool" in
        "last")
            # LAST is memory intensive
            echo $((available_memory / 4))  # 4GB per thread
            ;;
        "snpeff")
            # SnpEff uses Java heap
            echo $((available_memory / 10))  # 10GB per thread
            ;;
        "python")
            # Python scripts
            echo $((available_memory / 2))   # 2GB per thread
            ;;
        *)
            echo 4  # Default
            ;;
    esac
}

# ============================================================================
# ERROR HANDLING
# ============================================================================

handle_memory_error() {
    local error_msg="$1"
    echo "MEMORY ERROR: $error_msg"
    echo "Current memory usage:"
    free -h
    echo "Consider:"
    echo "  1. Reducing parallel processes"
    echo "  2. Processing data in smaller chunks"
    echo "  3. Increasing swap space"
    echo "  4. Using a machine with more memory"
    exit 1
}

# ============================================================================
# USAGE
# ============================================================================

usage() {
    echo "Memory Configuration for Phase 3 Annotation Pipeline"
    echo ""
    echo "Usage: source memory_config.sh"
    echo ""
    echo "Functions available:"
    echo "  check_memory_usage     - Check current memory usage"
    echo "  monitor_memory [sec]   - Monitor memory continuously"
    echo "  get_optimal_threads    - Get optimal thread count for tool"
    echo "  process_vcf_chunked   - Process VCF in chunks"
    echo "  cleanup_chunks        - Clean up chunk files"
    echo ""
    echo "Environment variables set:"
    echo "  JAVA_OPTS             - Java memory settings for SnpEff"
    echo "  MAX_MEMORY_GB         - Maximum memory limit"
}

# Show usage if script is run directly
if [ "${BASH_SOURCE[0]}" == "${0}" ]; then
    usage
fi
