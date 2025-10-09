# one may need to enable
#   Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
# to run this script
# or run: powershell -ExecutionPolicy Bypass -File .\run_tests_windows.ps1

# Clear the screen
Clear-Host

# Set environment variable
$env:OMP_NUM_THREADS = 1
$PY_COMMAND = "python"  # change to "python3" if needed

# Ensure log directory exists
$logDir = ".\ztest_logs"
if (-not (Test-Path $logDir)) {
    New-Item -ItemType Directory -Path $logDir | Out-Null
}

# Get a timestamp for log filename
$timestamp = Get-Date -Format "yyyy_MM_dd"

# Determine if script has arguments
if ($args.Count -eq 0) {
    Write-Host "Run all tests"
    $output = Join-Path $logDir ("$env:USERNAME" + "_" + "$env:COMPUTERNAME" + "_" + "$timestamp.log")
    Write-Host "Logging to $output"

    # Run the tests and log output
    & $PY_COMMAND "run_tests.py" $PY_COMMAND 2>&1 | Tee-Object -FilePath $output
} else {
    Write-Host "Run tests with arguments: $args"
    $output = Join-Path $logDir ("$env:USERNAME" + "_" + "$env:COMPUTERNAME" + "_" + "$timestamp-$($args[0]).log")
    Write-Host "Logging to $output"

    # Run the tests with arguments and log output
    & $PY_COMMAND "run_tests.py" $PY_COMMAND @args 2>&1 | Tee-Object -FilePath $output
}
