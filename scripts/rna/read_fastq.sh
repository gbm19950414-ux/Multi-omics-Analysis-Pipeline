#!/usr/bin/env bash
# 统一由 STAR 调用本脚本来解压 FASTQ（避免 PATH/引号/数组分词问题）
exec /usr/bin/gzip -dc "$@"
