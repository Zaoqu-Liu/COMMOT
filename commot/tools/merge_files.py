"""
合并优化代码到原文件，保持API完全兼容
"""

# 读取BACKUP（原文件）
with open('_spatial_communication_BACKUP.py', 'r') as f:
    original = f.read()

# 读取当前文件（已经是优化版本）
with open('_spatial_communication.py', 'r') as f:
    current = f.read()

# 策略：使用BACKUP作为基础，只替换关键部分
# 1. 替换CellCommunication类
# 2. 替换spatial_communication函数
# 3. 保留所有其他函数

# 从BACKUP提取所有函数定义
import re

# 找到所有需要保留的函数
functions_to_keep = [
    'communication_direction',
    'cluster_communication', 
    'cluster_communication_spatial_permutation',
    'cluster_position'
]

# 简单策略：使用BACKUP，然后插入优化的类和函数
print("正在合并文件...")

# 读取优化版本的关键部分
with open('_spatial_communication.py', 'r') as f:
    lines = f.readlines()

# 找到CellCommunicationHeavyOpt类的定义
opt_class_start = None
opt_class_end = None
for i, line in enumerate(lines):
    if 'class CellCommunicationHeavyOpt' in line:
        opt_class_start = i
    if opt_class_start and line.startswith('def ') and 'class' not in line:
        opt_class_end = i
        break

# 找到compute_sparse_distance_matrix
compute_func_start = None
compute_func_end = None
for i, line in enumerate(lines):
    if 'def compute_sparse_distance_matrix' in line:
        compute_func_start = i
    if compute_func_start and line.startswith('def ') and i > compute_func_start:
        compute_func_end = i
        break

# 找到spatial_communication函数
spatial_func_start = None
spatial_func_end = None
for i, line in enumerate(lines):
    if 'def spatial_communication(' in line and 'ultimate' not in line:
        spatial_func_start = i
    if spatial_func_start and i > spatial_func_start and (line.startswith('def ') or i == len(lines)-1):
        spatial_func_end = i if line.startswith('def ') else len(lines)
        break

print(f"找到优化类: {opt_class_start}-{opt_class_end}")
print(f"找到compute函数: {compute_func_start}-{compute_func_end}")
print(f"找到spatial函数: {spatial_func_start}-{spatial_func_end}")

# 现在构建最终文件
# 策略：从BACKUP开始，替换CellCommunication和spatial_communication

with open('_spatial_communication_BACKUP.py', 'r') as f:
    backup_lines = f.readlines()

# 找到BACKUP中的CellCommunication类
backup_class_start = None
backup_class_end = None
for i, line in enumerate(backup_lines):
    if 'class CellCommunication' in line:
        backup_class_start = i
    if backup_class_start and i > backup_class_start and line.startswith('    def run_cot_signaling'):
        # 找到run_cot_signaling的结束
        for j in range(i, len(backup_lines)):
            if backup_lines[j].startswith('def ') and 'class' not in backup_lines[j]:
                backup_class_end = j
                break
        break

print(f"BACKUP类位置: {backup_class_start}-{backup_class_end}")

# 构建最终文件
final_lines = []

# 1. 保留头部（imports等）
final_lines.extend(backup_lines[:backup_class_start])

# 2. 插入优化的类
if opt_class_start and opt_class_end:
    final_lines.extend(lines[opt_class_start:opt_class_end])
else:
    print("警告：未找到优化类，使用原类")
    final_lines.extend(backup_lines[backup_class_start:backup_class_end])

# 3. 插入compute_sparse_distance_matrix
if compute_func_start and compute_func_end:
    final_lines.extend(lines[compute_func_start:compute_func_end])

# 4. 找到BACKUP中的spatial_communication并替换
backup_spatial_start = None
backup_spatial_end = None
for i, line in enumerate(backup_lines):
    if 'def spatial_communication(' in line:
        backup_spatial_start = i
    if backup_spatial_start and i > backup_spatial_start:
        if (line.startswith('def ') and 'spatial_communication' not in line) or i == len(backup_lines)-1:
            backup_spatial_end = i if line.startswith('def ') else len(backup_lines)
            break

print(f"BACKUP spatial函数: {backup_spatial_start}-{backup_spatial_end}")

# 插入优化的spatial_communication
if spatial_func_start and spatial_func_end:
    final_lines.extend(lines[spatial_func_start:spatial_func_end])
else:
    print("警告：未找到优化函数，使用原函数")
    if backup_spatial_start and backup_spatial_end:
        final_lines.extend(backup_lines[backup_spatial_start:backup_spatial_end])

# 5. 保留BACKUP中spatial_communication之后的所有函数
if backup_spatial_end:
    final_lines.extend(backup_lines[backup_spatial_end:])

# 写入最终文件
with open('_spatial_communication.py', 'w') as f:
    f.writelines(final_lines)

print("✅ 文件合并完成！")
print(f"   总行数: {len(final_lines)}")
