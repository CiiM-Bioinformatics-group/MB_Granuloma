// 等待文档加载完成后执行
document.addEventListener('DOMContentLoaded', function() {
  var knowMoreBtn = document.getElementById('know-more-btn');
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

// **修正 `redirectToPage` 以正确处理 HTML 元素**
function redirectToPage(elementOrUrl) {
    if (typeof elementOrUrl === "string") {
        window.location.href = elementOrUrl;  // 直接跳转
    } else if (elementOrUrl instanceof HTMLElement) {
        var targetUrl = elementOrUrl.getAttribute("data-url");
        if (targetUrl) {
            window.location.href = targetUrl;
        } else {
            console.error("Missing data-url attribute on element:", elementOrUrl);
        }
    } else {
        console.error("Invalid parameter passed to redirectToPage:", elementOrUrl);
    }
}


// 从 cookie 中提取 CSRF token（Django 的默认机制）
function getCookie(name) {
  let cookieValue = null;
  if (document.cookie && document.cookie !== '') {
      const cookies = document.cookie.split(';');
      for (let i = 0; i < cookies.length; i++) {
          const cookie = cookies[i].trim();
          if (cookie.substring(0, name.length + 1) === (name + '=')) {
              cookieValue = decodeURIComponent(cookie.substring(name.length + 1));
              break;
          }
      }
  }
  return cookieValue;
}




  //添加contact

// 监听按钮点击事件
document.getElementById('sendButton').addEventListener('click', function() {
  // 获取表单元素
  const form = document.getElementById('contactForm');
  // 手动触发表单提交
  form.dispatchEvent(new Event('submit', { 'bubbles': true, 'cancelable': true }));
});

// 监听表单提交事件
document.getElementById('contactForm').addEventListener('submit', function(event) {
  event.preventDefault(); // 防止表单默认提交

  // 获取表单数据
  const companyName = document.getElementById('companyName').value;
  const fullName = document.getElementById('fullName').value;
  const emailAddress = document.getElementById('emailAddress').value;

  // 使用 fetch API 发送 POST 请求到后端 API
  // fetch('http://127.0.0.1:8000/api/contact/', {  // 替换为你的后端 API URL
  //     method: 'POST',
  //     headers: {
  //         'Content-Type': 'application/json',
  //     },
  //     body: JSON.stringify({
  //         company_name: companyName,
  //         full_name: fullName,
  //         email_address: emailAddress,
  //     }),
  // })

  const csrftoken = getCookie('csrftoken');

  fetch('http://127.0.0.1:8000/api/contact/', {
      method: 'POST',
      headers: {
          'Content-Type': 'application/json',
          'X-CSRFToken': csrftoken,  // 这是关键！
      },
      body: JSON.stringify({
          company_name: companyName,
          full_name: fullName,
          email_address: emailAddress,
      }),
  })



  .then(response => response.json())
  .then(data => {
      console.log('Success:', data);
      alert('Contact information submitted successfully!');
      // 你可以在这里清空表单或者重定向到另一个页面
      document.getElementById('contactForm').reset();
  })
  .catch((error) => {
      console.error('Error:', error);
      alert('Failed to submit contact information.');
  });
});

// 获取当前年份
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("footer-year").textContent = new Date().getFullYear();
});


