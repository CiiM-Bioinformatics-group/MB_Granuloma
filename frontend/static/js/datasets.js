
// // 等待文档加载完成后执行
// document.addEventListener('DOMContentLoaded', function() {
//   // 获取按钮元素
//   var knowMoreBtn = document.getElementById('Datasets-btn');

//   // 添加点击事件监听器
//   knowMoreBtn.addEventListener('click', function() {
//       // 跳转到 about.html 页面
//       window.location.href = 'lungs.html';
//   });
// });

// // 等待文档加载完成后执行
// document.addEventListener('DOMContentLoaded', function() {
//   // 获取按钮元素
//   var knowMoreBtn = document.getElementById('Sections-btn');

//   // 添加点击事件监听器
//   knowMoreBtn.addEventListener('click', function() {
//       // 跳转到 about.html 页面
//       window.location.href = 'lungs.html';
//   });
// });

// 等待文档加载完成后执行
document.addEventListener('DOMContentLoaded', function() {
  var knowMoreBtn = document.getElementById('Datasets-btn');
  if (knowMoreBtn) {
      knowMoreBtn.addEventListener('click', function() {
          // 获取 Django 传递的 URL
          var targetUrl = knowMoreBtn.getAttribute('data-url');
          if (targetUrl) {
              window.location.href = targetUrl;  // 直接跳转
          } else {
              console.error("Missing data-url attribute on button");
          }
      });
  }
});

// 等待文档加载完成后执行
document.addEventListener('DOMContentLoaded', function() {
  var knowMoreBtn = document.getElementById('Sections-btn');
  if (knowMoreBtn) {
      knowMoreBtn.addEventListener('click', function() {
          // 获取 Django 传递的 URL
          var targetUrl = knowMoreBtn.getAttribute('data-url');
          if (targetUrl) {
              window.location.href = targetUrl;  // 直接跳转
          } else {
              console.error("Missing data-url attribute on button");
          }
      });
  }
});


const header = document.querySelector("header");

    window.addEventListener("scroll", () => {
      const y = window.scrollY;

      if (y > 0 && y < 300) {
        header.classList.add("bg1");
      }

      if (y >= 300) {
        header.classList.add("bg2");
      }

      if (y == 0) {
        header.classList.remove("bg1", "bg2");
      }
    });

  
  const currentPath = window.location.pathname;
  document.querySelectorAll('.label-box').forEach(link => {
    if (currentPath.startsWith(link.getAttribute('href'))) {
      link.classList.add('active');
    }
  });

  document.addEventListener("DOMContentLoaded", function () {
    lucide.createIcons();
  });
  