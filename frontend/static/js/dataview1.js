// 自动识别 SCRIPT_NAME 路径前缀（如 "/mb-granuloma"）
// const SCRIPT_NAME = window.location.pathname.includes('/mb-granuloma') ? '/mb-granuloma' : '';

// 自动识别 SCRIPT_NAME 路径前缀（如 "/mb-granuloma/dataview1" 或 "/mb-granuloma/ACPMB003"）
const pathParts = window.location.pathname.split("/").filter(Boolean);
const SCRIPT_NAME = "/" + pathParts.slice(0, 2).join("/");
// 页面滚动联动目录高亮逻辑
const links = document.querySelectorAll("#list a[href^='#']");
const titles = [];

const highlight = id => {
    for (const link of links) {
        link.classList.remove("active");
    }
    id.classList.add("active");
};

for (const link of links) {
    link.addEventListener("click", () => highlight(link));
    const url = new URL(link.href);
    const title = document.querySelector(url.hash);
    titles.push(title);
}

const de = (func, dely = 100) => {
    let timer = null;
    return function () {
        clearTimeout(timer);
        timer = setTimeout(() => func.apply(this, arguments), dely);
    };
};

const scoll = de(() => {
    const range = [0, 300];
    const rects = titles.map(item => item.getBoundingClientRect());
    for (let i = 0; i < titles.length; i++) {
        const rect = rects[i];
        if (rect.top >= range[0] && rect.top <= range[1]) {
            highlight(links[i]);
            break;
        }
    }
});

window.addEventListener("scroll", scoll);

// 主逻辑
document.addEventListener("DOMContentLoaded", function () {
    const geneInput = document.getElementById("gene-selector");
    const suggestionBox = document.getElementById("gene-suggestions");
    const imgDisplay = document.getElementById("image-display");

    if (!imgDisplay) {
        console.warn("image-display 元素未找到！");
    }

    let allGenes = [];

    // 1️ 加载基因列表
    function loadGeneList() {
        fetch(`${SCRIPT_NAME}/api/gene_list/`)
            .then(response => response.json())
            .then(data => {
                if (data.genes) {
                    allGenes = data.genes;
                }
            })
            .catch(error => console.error("Error loading gene list:", error));
    }

    // 2️ 输入框联想
    geneInput.addEventListener("input", function () {
        const value = geneInput.value.toLowerCase();
        suggestionBox.innerHTML = "";

        if (!value) {
            suggestionBox.classList.add("hidden");
            return;
        }

        const filtered = allGenes.filter(g => g.toLowerCase().includes(value)).slice(0, 20);

        if (filtered.length === 0) {
            suggestionBox.classList.add("hidden");
            return;
        }

        filtered.forEach(gene => {
            const li = document.createElement("li");
            li.textContent = gene;
            li.addEventListener("click", function () {
                geneInput.value = gene;
                suggestionBox.classList.add("hidden");
                triggerGenePlot(gene);
            });
            suggestionBox.appendChild(li);
        });

        suggestionBox.classList.remove("hidden");
    });

    // 3️ 回车键绘图
    geneInput.addEventListener("keydown", function (e) {
        if (e.key === "Enter") {
            e.preventDefault();
            suggestionBox.classList.add("hidden");
            const selectedGene = geneInput.value;
            if (selectedGene) {
                triggerGenePlot(selectedGene);
            }
        }
    });

    // 4️ 绘图函数
    function triggerGenePlot(gene) {
        fetch(`${SCRIPT_NAME}/api/plot_gene_image/`, {
            method: "POST",
            headers: {
                "Content-Type": "application/json",
            },
            body: JSON.stringify({ gene })
        })
        .then(response => response.json())
        .then(data => {
            console.log("图像路径：", data.image_url);
            if (data.image_url && imgDisplay) {
                imgDisplay.src = data.image_url + "?t=" + new Date().getTime(); // 防止缓存
            } else {
                console.error("Failed to load image:", data.error);
            }
        });
    }

    // 5️ 页面加载时初始化
    loadGeneList();
});