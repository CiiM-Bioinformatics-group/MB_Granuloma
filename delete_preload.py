import os


apps = [
    'ACPMB002', 'ACPMB003', 'ACPMB004', 'ACPMB005', 'ACPMB006', 'ACPMB007',
    'ACPMB008', 'ACPMB009', 'ACPMB010', 'ACPMB011', 'ACPMB013', 'ACPMB014', 'ACPMB016',
    'NSG005', 'NSG013', 'NSG014', 'NSG017', 'NSG019', 'NSG020_1', 'NSG020_2', 'NSG020_3',
    'NSG023', 'NSG024', 'YTB001', 'YTB005', 'YTB006', 'YTB008', 'YTB009', 'YTB010',
    'YTB012', 'YTB015', 'YTB019'
]

# add dataview1 but skip it
apps_with_dataview1 = apps + ["dataview1"]

for app in apps_with_dataview1:
    preload_path = os.path.join(app, "preload.py")
    if app == "dataview1":
        print(f"🔒 Skipping: {preload_path}")
        continue
    if os.path.exists(preload_path):
        os.remove(preload_path)
        print(f"🗑 Deleted: {preload_path}")
    else:
        print(f"❌ Not found: {preload_path}")
