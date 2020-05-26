<template>
  <div>
      <div class="LoginTop">
          <div>1</div>
          <div>{{middTop}}</div>
          <div>1</div>
      </div>
  </div>
</template>

<script>
export default {
    props:['middTop']
}
</script>

<style>  */ less存在问题
.LoginTop{
    height: 50px;
    background-color: white;
    display: flex;
    /* div{
        flex: 1;
        display: flex;
        justify-content: center;
        align-items: center;
    }; */
}

</style>