

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

//     // 等待文档加载完成后执行
// document.addEventListener('DOMContentLoaded', function() {
//   // 获取按钮元素
//   var knowMoreBtn = document.getElementById('Datasets-bth');

//   // 添加点击事件监听器
//   knowMoreBtn.addEventListener('click', function() {
//       // 跳转到 about.html 页面
//       window.location.href = 'lungs.html';
//   });

//   var knowMoreBtn = document.getElementById('Sectionss-bth');

//   // 添加点击事件监听器
//   knowMoreBtn.addEventListener('click', function() {
//       // 跳转到 about.html 页面
//       window.location.href = 'lungs.html';
//   });
// });

// // 等待文档加载完成后执行
// document.addEventListener('DOMContentLoaded', function() {
//   var knowMoreBtn = document.getElementById('know-more-btn');
//   if (knowMoreBtn) {
//       knowMoreBtn.addEventListener('click', function() {
//           // 获取 Django 传递的 URL
//           var targetUrl = knowMoreBtn.getAttribute('data-url');
//           if (targetUrl) {
//               window.location.href = targetUrl;  // 直接跳转
//           } else {
//               console.error("Missing data-url attribute on button");
//           }
//       });
//   }
// });