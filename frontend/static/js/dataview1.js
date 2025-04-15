
const SCRIPT_NAME = '/mb-granuloma';

      const links = document.querySelectorAll("#list a[href^='#']");
  
      const titles = [];
  
      const highlight = id => {
        for (const link of links) {
          link.classList.remove("active");
        }
  
        id.classList.add("active");
      };
  
      for (const link of links) {
        link.addEventListener("click", () => {
          highlight(link);
        });
  
        const url = new URL(link.href);
  
        const title = document.querySelector(url.hash);
  
        titles.push(title);
      }
  
      const de = (func, dely = 100) => {
        let timer = null;
  
        return function () {
          clearTimeout(timer);
          timer = setTimeout(() => {
            func.apply(this, arguments);
          }, dely);
        };
      };
  
      const scoll = de(() => {
        const range = [0, 300];
        const rects = titles.map(item => item.getBoundingClientRect());
  
        for (let i = 0; i < titles.length; i++) {
          const title = titles[i];
          const rect = rects[i];
  
          if (rect.top >= range[0] && rect.top <= range[1]) {
            highlight(links[i]);
            break;
          }
        }
      });
  
      window.addEventListener("scroll", scoll);

    //   // 以下是为了可视化visium数据
    //   document.addEventListener("DOMContentLoaded", function () {
    //     const geneSelect = document.querySelector("#gene-selector");
    //     const imgDisplay = document.getElementById("image-display");
    
    //     // 1. 获取基因列表

    //       function loadGeneList() {
    //         fetch("/dataview1/api/gene_list/")
    //         .then(response => response.json())
    //         .then(data => {
    //             if (data.genes) {
    //                 const dataList = document.getElementById("gene-list");
    //                 dataList.innerHTML = ""; // 清空旧的 option
    //                 data.genes.forEach(gene => {
    //                     const option = document.createElement("option");
    //                     option.value = gene;
    //                     dataList.appendChild(option);
    //                 });
    //             }
    //         })
    //         .catch(error => console.error("Error loading gene list:", error));
    //     }
      
    
    //     // 2. 用户选择基因后，调用后端生成图像
    //     geneSelect.addEventListener("change", function () {
    //         const selectedGene = geneSelect.value;
    //         if (selectedGene) {
    //             fetch("/dataview1/api/plot_gene_image/", {
    //                 method: "POST",
    //                 headers: {
    //                     "Content-Type": "application/json",
    //                 },
    //                 body: JSON.stringify({ gene: selectedGene })
    //             })
    //             .then(response => response.json())
    //             .then(data => {
    //               console.log("图像路径：", data.image_url);  //  新加的行
    //                 if (data.image_url) {
    //                     imgDisplay.src = data.image_url + "?t=" + new Date().getTime();  // 加时间戳防缓存
    //                 } else {
    //                     console.error("Failed to load image:", data.error);
    //                 }
    //             });
    //         }
    //     });
        
    
    //     // 3. 页面初始化加载基因列表
    //     loadGeneList();
    // });


  //   document.addEventListener("DOMContentLoaded", function () {
  //     const geneInput = document.getElementById("gene-selector");
  //     const imgDisplay = document.getElementById("image-display");
  
  //     function loadGeneList() {
  //         fetch("/dataview1/api/gene_list/")
  //         .then(response => response.json())
  //         .then(data => {
  //             if (data.genes) {
  //                 const dataList = document.getElementById("gene-list");
  //                 dataList.innerHTML = "";
  //                 data.genes.forEach(gene => {
  //                     const option = document.createElement("option");
  //                     option.value = gene;
  //                     dataList.appendChild(option);
  //                 });
  //             }
  //         })
  //         .catch(error => console.error("Error loading gene list:", error));
  //     }
  
  //     geneInput.addEventListener("change", function () {
  //         const selectedGene = geneInput.value;
  //         if (selectedGene) {
  //             fetch("/dataview1/api/plot_gene_image/", {
  //                 method: "POST",
  //                 headers: {
  //                     "Content-Type": "application/json",
  //                 },
  //                 body: JSON.stringify({ gene: selectedGene })
  //             })
  //             .then(response => response.json())
  //             .then(data => {
  //                 console.log("图像路径：", data.image_url);
  //                 if (data.image_url) {
  //                     imgDisplay.src = data.image_url + "?t=" + new Date().getTime();
  //                 } else {
  //                     console.error("Failed to load image:", data.error);
  //                 }
  //             });
  //         }
  //     });
  
  //     loadGeneList();
  // });

  document.addEventListener("DOMContentLoaded", function () {
    const geneInput = document.getElementById("gene-selector");
    const suggestionBox = document.getElementById("gene-suggestions");
    const imgDisplay = document.getElementById("image-display");
    let allGenes = [];

    // 1️ 载入 gene 列表（放入 JS 变量中）
    function loadGeneList() {
      fetch(`${SCRIPT_NAME}/dataview1/api/gene_list/`)
            .then(response => response.json())
            .then(data => {
                if (data.genes) {
                    allGenes = data.genes;
                }
            })
            .catch(error => console.error("Error loading gene list:", error));
    }

    // 2️ 输入框监听用户输入，展示匹配建议
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

    // 3️ 回车时直接触发绘图
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

    // 4️ 核心函数：发送请求并更新图片
    function triggerGenePlot(gene) {
      fetch(`${SCRIPT_NAME}/dataview1/api/gene_list/`, {
            method: "POST",
            headers: {
                "Content-Type": "application/json",
            },
            body: JSON.stringify({ gene })
        })
        .then(response => response.json())
        .then(data => {
            console.log("图像路径：", data.image_url);
            if (data.image_url) {
                imgDisplay.src = data.image_url + "?t=" + new Date().getTime(); // 防止缓存
            } else {
                console.error("Failed to load image:", data.error);
            }
        });
    }

    // 5️ 页面加载时初始化
    loadGeneList();
});

  
    